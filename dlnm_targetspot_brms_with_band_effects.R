# DLNM + Bayesian Beta regression (brms) for Target Spot final severity (soybean)
# -------------------------------------------------------------------------
# Inputs (created from ts_check.xlsx and weather_nasa.xlsx):
#   - targetspot_epidemics_clean.csv
#   - targetspot_weather_long_0_90dap.csv
#   - targetspot_assessment_endpoints_80_85_90dap.csv
#
# This script:
#   1) builds DLNM cross-basis design matrices for selected weather variables
#   2) fits a hierarchical Beta regression in brms with site-year random effects
#   3) summarizes lag-band cumulative effects (0–40 vs 41–LAG_MAX) from the posterior
#
# Key references:
# - Gasparrini A (2011) Distributed lag linear and non-linear models in R: the package dlnm.
#   Journal of Statistical Software 43(8):1–20.
# - Ferrari SLP, Cribari-Neto F (2004) Beta regression for modelling rates and proportions.
#   Journal of Applied Statistics 31:799–815.
# - Gelman A, Hill J (2007) Data Analysis Using Regression and Multilevel/Hierarchical Models. CUP.
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(dlnm)
  library(brms)
})

# ------------------------------- PARAMETERS ----------------------------------

END_DAP   <- 85      # choose: 80, 85, or 90 (endpoint = planting + END_DAP days)
LAG_MAX   <- 70      # max lag (days before END_DAP) included in DLNM
DF_VAR    <- 4       # exposure-response complexity (natural splines)
DF_LAG    <- 4       # lag-response complexity (natural splines)
SEPARATOR <- LAG_MAX # number of NA separators between epidemics when building template cb

# Weather variables to include (must exist in weather_long csv):
# - tmean is computed in the preprocessing: (T2M_MIN + T2M_MAX)/2 if available, else T2M
# - rh: RH2M
# - rain: PRECTOTCORR
VARS <- c("tmean", "vpd")

# lag-band definitions (in days before END_DAP)
LAG_BAND_A <- 0:40
LAG_BAND_B <- 41:LAG_MAX

# posterior sampling for effect summaries
NDRAWS <- 600   # number of posterior draws used for lag-band summaries (thinned)
AT_PROBS <- c(0.50, 0.75)  # contrasts: 75th vs 50th percentile ("high vs typical")

# brms settings
CHAINS <- 4
ITER   <- 4000
CORES  <- 4

# ------------------------------- READ DATA -----------------------------------

epi  <- read.csv("data/targetspot_epidemics_clean.csv", stringsAsFactors = FALSE)
wx   <- read.csv("data/targetspot_weather_long_0_90dap.csv", stringsAsFactors = FALSE)
ends <- read.csv("data/targetspot_assessment_endpoints_80_85_90dap.csv", stringsAsFactors = FALSE)

stopifnot(END_DAP %in% unique(ends$end_dap))

# endpoint date per epidemic
end_use <- ends %>% filter(end_dap == END_DAP) %>% dplyr::select(epi_id, assessment_date_assumed)

# analysis set: non-missing severity AND has weather
epi_use <- epi %>%
  inner_join(end_use, by = "epi_id") %>%
  filter(!is.na(sev)) %>%
  mutate(siteyear = ifelse(is.na(siteyear) | siteyear == "", paste0(location, "_", year), siteyear))

# keep weather within 0..END_DAP, and only for epidemics in epi_use
wx_use <- wx %>%
  filter(dpp >= 0, dpp <= END_DAP) %>%
  semi_join(epi_use, by = "epi_id") %>%
  arrange(epi_id, dpp)


es <- function(T) {
  0.6108 * exp((17.27 * T) / (T + 237.3))
}

# Ensure VPD exists and is consistent with tmean + rh
wx_use <- wx_use %>%
  mutate(
    vpd = es(tmean) * (1 - rh / 100)
  )


# completeness: need all days 0..END_DAP per epidemic for clean lags
comp <- wx_use %>%
  group_by(epi_id) %>%
  summarise(n_days = n_distinct(dpp), .groups="drop") %>%
  mutate(ok = n_days == (END_DAP + 1))

message("Epidemics with complete daily weather (0..END_DAP): ", sum(comp$ok), " / ", nrow(comp))

epi_use <- epi_use %>% semi_join(comp %>% filter(ok), by="epi_id")
wx_use  <- wx_use  %>% semi_join(comp %>% filter(ok), by="epi_id")

# ------------------------ BUILD TEMPLATE CROSSBASIS ---------------------------
# We build one template crossbasis per variable using pooled data separated by NA blocks.
# This prevents lag structures from "bleeding" across epidemics, while letting dlnm pick
# consistent spline bases for exposure and lag.

build_pooled_series <- function(wx_long, var, sep_n = SEPARATOR) {
  # returns a numeric vector with NA separators between epidemics
  ids <- unique(wx_long$epi_id)
  out <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    v <- wx_long %>% filter(epi_id == ids[i]) %>% arrange(dpp) %>% pull(.data[[var]])
    out[[i]] <- c(v, rep(NA_real_, sep_n))
  }
  unlist(out)
}

# Create template crossbasis objects (one per variable)
cb_templates <- list()
for (v in VARS) {
  if (!v %in% names(wx_use)) stop("Weather variable not found: ", v)
  x_pool <- build_pooled_series(wx_use, v, sep_n = SEPARATOR)
  cb_templates[[v]] <- crossbasis(x_pool, lag = LAG_MAX,
                                 argvar = list(fun = "ns", df = DF_VAR),
                                 arglag = list(fun = "ns", df = DF_LAG))
}

# ---------------------- BUILD DESIGN MATRIX (ONE ROW / EPIDEMIC) -------------
# For each epidemic, we compute a crossbasis on its own time series and extract
# the LAST ROW (at END_DAP), which summarizes exposures over lags 0..LAG_MAX.

extract_last_cb_row <- function(x, cb_template) {
  cb <- crossbasis(x, lag = LAG_MAX,
                   argvar = attr(cb_template, "argvar"),
                   arglag = attr(cb_template, "arglag"))
  as.numeric(cb[length(x), ])
}

build_design_for_var <- function(wx_long, var, cb_template, prefix) {
  X <- wx_long %>%
    group_by(epi_id) %>%
    summarise(cb = list(extract_last_cb_row(.data[[var]], cb_template)),
              .groups = "drop")
  # name columns deterministically
  p <- length(X$cb[[1]])
  nm <- paste0(prefix, seq_len(p))
  X %>% mutate(cb = lapply(cb, setNames, nm)) %>% unnest_wider(cb)
}

X_list <- list()
for (v in VARS) {
  X_list[[v]] <- build_design_for_var(wx_use, v, cb_templates[[v]], prefix = paste0("cb_", v, "_"))
}
library(tidyverse)
X <- reduce(X_list, left_join, by = "epi_id")

dat <- epi_use %>%
  select(epi_id, sev, siteyear, year, location, state, lat, lon, planting_date, assessment_date_assumed) %>%
  left_join(X, by = "epi_id")

# ------------------------------- FIT MODEL -----------------------------------
# Beta regression for y in (0,1). Your sev was pre-adjusted into (0,1) already.
# Random intercept for siteyear (site-year environment).

# Build formula: include all cb columns
cb_terms <- paste(names(dat)[grepl("^cb_", names(dat))], collapse = " + ")
fml <- as.formula(paste0("sev ~ 1 + (1|siteyear) + ", cb_terms))

priors <- c(
  prior(normal(0, 0.7), class = "b"),        # a bit tighter than N(0,1)
  prior(exponential(1), class = "phi"),
  prior(student_t(3, 0, 2.5), class = "sd")
)



fit <- brm(
  formula = fml,
  data    = dat,
  family  = Beta(),
  prior   = priors,
  chains  = CHAINS,
  iter    = ITER,
  cores   = CORES,
  seed    = 20251218,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

pp_check(fit, type="dens_overlay")
bayes_R2(fit)
predict(fit)
print(fit)
summary(fit)

# ----------------------- LAG-BAND EFFECT SUMMARIES ---------------------------
# Goal: summarize "cumulative" effects in two lag bands:
#   A: lags 0–40 days before END_DAP (secondary cycles / epidemic build-up)
#   B: lags 41–LAG_MAX (initiation / early infections)
#
# We compute, for each variable separately, the cumulative exposure effect in each lag band
# for a contrast between AT_PROBS[2] vs AT_PROBS[1] (e.g., 75th vs 50th percentile).
#
# Approach:
#   - Use the TEMPLATE crossbasis (pooled) to build a crosspred object for each posterior draw
#     using coef=... for the relevant cb columns.
#   - crosspred returns matfit (at x lag) on the linear predictor scale.
#   - Summation across lag columns within a band yields band-specific cumulative effect.
#
# NOTE: This is an effect on the linear predictor (logit mean). We report:
#   - delta_eta_band (high - typical) and its exp() as an odds-ratio on the mean scale.

# helper: get posterior draws for coefficients of a given variable
get_draws_for_prefix <- function(fit, prefix) {
  draws <- as_draws_df(fit)
  cols <- grep(paste0("^b_", prefix), names(draws), value = TRUE)
  if (length(cols) == 0) stop("No coefficients found for prefix: ", prefix)
  draws %>% select(all_of(cols))
}

# helper: compute band effects for one variable
band_effects_one_var <- function(var, cb_template, fit, at_probs = AT_PROBS,
                                 lagA = LAG_BAND_A, lagB = LAG_BAND_B, ndraws = NDRAWS) {

  # exposure grid for crosspred
  x_pool <- attr(cb_template, "range") # might be NULL; safer to use data-based quantiles
  # derive at values from wx_use pooled distribution
  
  
  x_all <- wx_use %>% pull(.data[[var]])
  at_vals <- as.numeric(quantile(x_all, probs = at_probs, na.rm = TRUE))
  names(at_vals) <- paste0("p", at_probs*100)

  # posterior draws (thin to ndraws)
  prefix <- paste0("cb_", var, "_")
  d <- get_draws_for_prefix(fit, prefix)
  if (nrow(d) > ndraws) {
    set.seed(20251218)
    idx <- sort(sample(seq_len(nrow(d)), ndraws))
    d <- d[idx, , drop = FALSE]
  }

  # map brms coefficient names to cb column order
  # brms names: b_cb_var_1, b_cb_var_2, ...
  cb_cols <- grep(paste0("^", prefix), colnames(dat), value = TRUE)
  # ensure deterministic order by numeric suffix
  ord <- order(as.integer(gsub(paste0("^", prefix), "", cb_cols)))
  cb_cols <- cb_cols[ord]
  bnames <- paste0("b_", cb_cols)

  # pre-allocate
  out <- vector("list", nrow(d))

  for (i in seq_len(nrow(d))) {
    co <- as.numeric(d[i, bnames])
    names(co) <- bnames

    cp <- crosspred(cb_template,
                    coef = co,
                    vcov = diag(0, nrow = length(co)), # point estimate per draw
                    at   = at_vals,
                    cen  = at_vals[1],
                    bylag = 1)

    # cp$matfit is (length(at) x (LAG_MAX+1)) typically, lag columns 0..LAG_MAX
    mf <- cp$matfit
    # ensure lag indices correspond to 0..LAG_MAX
    lag_names <- colnames(mf)
    # Some versions use "lag0" etc; handle both
    lag_idx <- suppressWarnings(as.integer(gsub("lag", "", lag_names)))
    if (anyNA(lag_idx)) lag_idx <- seq(0, ncol(mf)-1)

    # locate high and typical rows
    r_typ <- 1
    r_hi  <- 2

    # cumulative per band for each at value are rowSums across lag subset
    colsA <- which(lag_idx %in% lagA)
    colsB <- which(lag_idx %in% lagB)

    # effect is already centered at at_vals[1]; so mf[r_hi, ] is delta_eta at each lag
    dA <- sum(mf[r_hi, colsA], na.rm = TRUE)
    dB <- sum(mf[r_hi, colsB], na.rm = TRUE)

    out[[i]] <- data.frame(
      var = var,
      end_dap = END_DAP,
      lag_max = LAG_MAX,
      at_typ = at_vals[1],
      at_hi  = at_vals[2],
      delta_eta_bandA = dA,
      delta_eta_bandB = dB,
      OR_bandA = exp(dA),
      OR_bandB = exp(dB)
    )
  }

  bind_rows(out)
}

band_draws <- bind_rows(lapply(VARS, function(v) band_effects_one_var(v, cb_templates[[v]], fit)))

# summarize posterior
band_summary <- band_draws %>%
  pivot_longer(cols = c(delta_eta_bandA, delta_eta_bandB, OR_bandA, OR_bandB),
               names_to = "metric", values_to = "value") %>%
  group_by(var, end_dap, lag_max, metric) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    lo95   = quantile(value, 0.025, na.rm = TRUE),
    hi95   = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(band_draws,   file = paste0("dlnm_band_effects_draws_END", END_DAP, "_L", LAG_MAX, ".csv"), row.names = FALSE)
write.csv(band_summary, file = paste0("dlnm_band_effects_summary_END", END_DAP, "_L", LAG_MAX, ".csv"), row.names = FALSE)

print(band_summary)

message("Saved band effect draws and summaries for END_DAP=", END_DAP, ", LAG_MAX=", LAG_MAX)

# ---------------------------- OPTIONAL: QUICK PLOT ----------------------------
# Simple forest-style plot for ORs
# (You can replace with ggplot2 if desired; kept dependency-light.)
# Example usage:
#   bs <- read.csv("dlnm_band_effects_summary_END85_L70.csv")
#   subset(bs, metric %in% c("OR_bandA","OR_bandB"))




library(ggplot2)
library(dplyr)

bs <- band_summary %>%
  filter(metric %in% c("OR_bandA", "OR_bandB")) %>%
  mutate(
    band = ifelse(metric == "OR_bandA",
                  "0–40 days before R6 (epidemic development)",
                  "41–70 days before R6 (initiation)"),
    var = factor(var, levels = c("tmean", "vpd"),
                 labels = c("Mean temperature", "VPD"))
  )

p <- ggplot(bs, aes(x = median, y = var)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95),
                 height = 0.15, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 4),
    limits = c(0.4, 6)
  ) +
  facet_wrap(~ band, ncol = 1) +
  labs(
    x = "Odds ratio (75th vs 50th percentile)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


print(p)







dat_dec <- dat %>%
  mutate(
    y20 = as.integer(sev > 0.35)   # 1 = severidade “alta” (aciona intensificação)
  )
table(dat_dec$y20)


get_lag_weights <- function(var, fit, cb_templates, wx_use, dat,
                            at_probs = c(0.50, 0.75), lagA = 0:40) {
  
  cb <- cb_templates[[var]]
  
  # contraste (p75 vs p50)
  x_all   <- wx_use[[var]]
  at_vals <- as.numeric(quantile(x_all, probs = at_probs, na.rm = TRUE))
  
  # coeficientes (mediana posterior) na ordem correta
  prefix <- paste0("cb_", var, "_")
  cb_cols <- grep(paste0("^", prefix), colnames(dat), value = TRUE)
  cb_cols <- cb_cols[order(as.integer(sub(prefix, "", cb_cols)))]
  bnames  <- paste0("b_", cb_cols)
  
  draws <- posterior::as_draws_df(fit)
  beta_med <- apply(draws[, bnames, drop = FALSE], 2, median)
  beta_med <- as.numeric(beta_med)
  
  stopifnot(length(beta_med) == ncol(cb))
  
  V0 <- matrix(0, nrow = length(beta_med), ncol = length(beta_med))
  
  cp <- dlnm::crosspred(
    cb, coef = beta_med, vcov = V0,
    at = at_vals, cen = at_vals[1], bylag = 1
  )
  
  # efeito por lag (não cumulativo) no contraste p75 vs p50:
  eff_lag <- as.numeric(cp$matfit[2, ])  # lags 0..LAG_MAX
  
  # manter apenas lags da banda A e normalizar (opcional)
  w <- eff_lag[lagA + 1]                # +1 pq R index
  w <- w / sum(abs(w))                  # escala comparável entre variáveis
  
  data.frame(lag = lagA, w = w)
}

w_t <- get_lag_weights("tmean", fit, cb_templates, wx_use, dat)
w_v <- get_lag_weights("vpd",   fit, cb_templates, wx_use, dat)


make_index <- function(wx_use, var, weights, END_DAP = 85) {
  # índice no momento END_DAP: usa os valores diários em (END_DAP - lag)
  # weights: data.frame(lag, w)
  
  wx_use %>%
    mutate(day = dpp) %>%
    inner_join(weights, by = character()) %>%  # truque: expand
    mutate(day_use = END_DAP - lag) %>%
    group_by(epi_id, lag) %>%
    summarise(x = .data[[var]][match(day_use, day)], w = first(w), .groups="drop") %>%
    group_by(epi_id) %>%
    summarise(idx = sum(w * x, na.rm = TRUE), .groups="drop")
}

idx_t <- make_index(wx_use, "tmean", w_t, END_DAP = END_DAP) %>% rename(idx_tmean = idx)
idx_v <- make_index(wx_use, "vpd",   w_v, END_DAP = END_DAP) %>% rename(idx_vpd   = idx)

dat_dec <- dat_dec %>%
  left_join(idx_t, by = "epi_id") %>%
  left_join(idx_v, by = "epi_id") %>%
  mutate(
    idx = scale(idx_tmean) + scale(idx_vpd)   # índice combinado (simples)
  )

library(brms)

fit20 <- brm(
  y20 ~ 1 + idx + (1 | siteyear),
  data = dat_dec,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(student_t(3, 0, 2.5), class = "sd")
  ),
  chains = 4, iter = 4000, cores = 4,
  seed = 20251218,
  control = list(adapt_delta = 0.95)
)

summary(fit20)

dat_dec$prob_high <- fitted(fit20, scale = "response")[, "Estimate"]




make_index_prospective <- function(wx_use, var, weights,
                                   TODAY_DAP,
                                   END_DAP = 85) {
  
  lag_max_eff <- min(max(weights$lag), END_DAP - TODAY_DAP)
  
  wx_use %>%
    filter(dpp <= TODAY_DAP) %>%   # <<< só clima observado até hoje
    mutate(day = dpp) %>%
    inner_join(weights %>% filter(lag <= lag_max_eff),
               by = character()) %>%
    mutate(day_use = TODAY_DAP - lag) %>%
    group_by(epi_id, lag) %>%
    summarise(
      x = .data[[var]][match(day_use, day)],
      w = first(w),
      .groups = "drop"
    ) %>%
    group_by(epi_id) %>%
    summarise(idx = sum(w * x, na.rm = TRUE), .groups = "drop")
}


TODAY_DAP <- 60

idx_t <- make_index_prospective(wx_use, "tmean", w_t,
                                TODAY_DAP = TODAY_DAP,
                                END_DAP = END_DAP) %>%
  rename(idx_tmean = idx)

idx_v <- make_index_prospective(wx_use, "vpd", w_v,
                                TODAY_DAP = TODAY_DAP,
                                END_DAP = END_DAP) %>%
  rename(idx_vpd = idx)

dat_dec <- dat %>%
  mutate(y20 = as.integer(sev > 0.35)) %>%
  left_join(idx_t, by = "epi_id") %>%
  left_join(idx_v, by = "epi_id") %>%
  mutate(idx = scale(idx_tmean) + scale(idx_vpd))


fit20 <- brm(
  y20 ~ 1 + idx + (1 | siteyear),
  data = dat_dec,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(student_t(3, 0, 2.5), class = "sd")
  ),
  chains = 4, iter = 4000, cores = 4,
  seed = 20251218,
  control = list(adapt_delta = 0.95)
)


library(dlnm)
library(posterior)
library(dplyr)
library(ggplot2)

var  <- "tmean"
cb   <- cb_templates[[var]]
lagA <- 0:40

# grade de temperatura (ex.: percentis 5 a 95)
t_all <- wx_use[[var]]
at_grid <- as.numeric(quantile(t_all, probs = seq(0.05, 0.95, by = 0.01), na.rm = TRUE))
at_grid <- sort(unique(at_grid))

# escolha do centro (referência): mediana
cen_val <- as.numeric(quantile(t_all, probs = 0.50, na.rm = TRUE))

# coeficientes do brms na ordem correta
prefix <- paste0("cb_", var, "_")
cb_cols <- grep(paste0("^", prefix), colnames(dat), value = TRUE)
cb_cols <- cb_cols[order(as.integer(sub(prefix, "", cb_cols)))]
bnames  <- paste0("b_", cb_cols)

draws <- as_draws_df(fit)

# usar NDRAWS draws (para IC)
NDRAWS <- 600
set.seed(1)
idx <- sample(seq_len(nrow(draws)), min(NDRAWS, nrow(draws)))
draws_sub <- draws[idx, bnames, drop = FALSE]

# função: para um draw, devolve RR (exp(delta_eta)) por temperatura na banda A
one_draw_rr <- function(beta_vec) {
  beta_vec <- as.numeric(beta_vec)
  V0 <- matrix(0, length(beta_vec), length(beta_vec))
  cp <- crosspred(
    cb,
    coef = beta_vec,
    vcov = V0,
    at   = at_grid,
    cen  = cen_val,
    bylag = 1
  )
  mf <- cp$matfit  # (length(at) x (LAG_MAX+1))
  lag_names <- colnames(mf)
  lag_idx <- suppressWarnings(as.integer(gsub("lag", "", lag_names)))
  if (anyNA(lag_idx)) lag_idx <- 0:(ncol(mf)-1)
  
  colsA <- which(lag_idx %in% lagA)
  
  # efeito cumulativo (logit scale) por T: soma nos lags da banda A
  etaA <- rowSums(mf[, colsA, drop = FALSE], na.rm = TRUE)
  
  # converte para razão de odds na escala do preditor médio (interpretável como “RR” relativo)
  exp(etaA)
}

RR_mat <- apply(draws_sub, 1, one_draw_rr)  # colunas = draws; linhas = at_grid

# sumarizar posterior por temperatura
rr_sum <- data.frame(
  tmean = at_grid,
  RR_med = apply(RR_mat, 1, median, na.rm = TRUE),
  RR_lo  = apply(RR_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
  RR_hi  = apply(RR_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
)

# faixa “mais favorável”: onde RR_med é máximo
t_best <- rr_sum$tmean[which.max(rr_sum$RR_med)]
t_best

ggplot(rr_sum, aes(x = tmean, y = RR_med)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_ribbon(aes(ymin = RR_lo, ymax = RR_hi), alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(
    x = "Temperatura média (°C)",
    y = "Razão relativa (75% vs referência, acumulada em lags 0–40)",
    title = "Faixa de temperatura mais favorável (efeito cumulativo na janela 0–40 dias antes de R6)"
  ) +
  theme_bw()



R.version.string
Sys.which("make")
pkgbuild::has_build_tools(debug = TRUE)



