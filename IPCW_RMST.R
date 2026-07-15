# IPCW RMST policy evaluation for the single-partition setting.
#
# This file implements a fixed-policy IPCW estimator of truncated mean survival
# time up to tau using the long-format data generated in SimData.R.

library(dplyr)
library(ranger)
library(survival)


# Validates required columns and basic coding for the IPCW estimator input.
validate_ipcw_data <- function(long_dat, history_vars) {
  required_cols <- unique(c("id", "k", "Ak", "Xk", "delta_k", "Bk", history_vars))
  missing_cols <- setdiff(required_cols, names(long_dat))
  if (length(missing_cols) > 0L) {
    stop("long_dat is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!all(long_dat$delta_k %in% c(0, 1, 0L, 1L, NA))) {
    stop("delta_k must be coded as 0/1.")
  }

  invisible(TRUE)
}

# Evaluates a right-continuous step survival curve at arbitrary times.
step_eval_surv_ipcw <- function(time_grid, surv_vec, t) {
  if (length(time_grid) == 0L) {
    return(rep(1, length(t)))
  }

  idx <- findInterval(t, time_grid)
  out <- ifelse(idx <= 0L, 1, surv_vec[pmax(idx, 1L)])
  out[t < 0] <- 1
  pmin(pmax(out, 0), 1)
}

# Computes the effective sample size of a nonnegative weight vector.
effective_sample_size_ipcw <- function(w) {
  w_pos <- w[is.finite(w) & w > 0]
  if (length(w_pos) == 0L) {
    return(0)
  }
  sum(w_pos)^2 / sum(w_pos^2)
}

# Adds row-level policy actions and policy-match indicators.
add_policy_match_ipcw <- function(dat, policy_fun) {
  dat$pi_hat <- as.integer(policy_fun(dat))
  dat$M_row <- as.integer(dat$Ak == dat$pi_hat)
  dat
}

# Fits a binary propensity model and returns observed-treatment probabilities.
estimate_propensity_ipcw <- function(dat, history_vars) {
  f <- as.formula(paste("Ak ~", paste(history_vars, collapse = " + ")))
  fit <- glm(f, data = dat, family = binomial())
  p1 <- as.numeric(predict(fit, type = "response"))
  p1 <- pmin(pmax(p1, 1e-4), 1 - 1e-4)
  p_obs <- ifelse(dat$Ak == 1L, p1, 1 - p1)
  list(fit = fit, p_obs = p_obs)
}

# Defines non-administrative censoring at the visit level before tau.
define_censoring_event_ipcw <- function(dat, tau) {
  dat$visit_stop_time <- pmin(dat$Bk + dat$Xk, tau)
  dat$cens_event_k <- as.integer(dat$delta_k == 0L & (dat$Bk + dat$Xk) < tau)
  dat
}

# Fits a censoring survival forest on the long-format visit data.
fit_censoring_forest_ipcw <- function(dat, history_vars,
                                      num.trees = 300, min.node.size = 10) {
  f <- as.formula(
    paste0("Surv(Xk, cens_event_k) ~ ", paste(c(history_vars, "Ak"), collapse = " + "))
  )
  ranger(
    formula = f,
    data = dat,
    num.trees = num.trees,
    min.node.size = min.node.size,
    respect.unordered.factors = "order",
    seed = 3030
  )
}

# Predicts row-level censoring survival probabilities at each observed visit length.
predict_censor_survival_ipcw <- function(fit, newdata) {
  pr <- predict(fit, data = newdata)
  base_grid <- pr$unique.death.times
  S <- pr$survival

  g_hat <- numeric(nrow(newdata))
  for (i in seq_len(nrow(newdata))) {
    t_eval <- max(newdata$Xk[i] - 1e-8, 0)
    g_hat[i] <- step_eval_surv_ipcw(base_grid, S[i, ], t_eval)
  }
  pmax(g_hat, 1e-4)
}

# Aggregates row-level quantities into subject-level IPCW contributions.
build_subject_level_ipcw <- function(dat, tau) {
  dat %>%
    arrange(id, k) %>%
    group_by(id) %>%
    summarise(
      n_visits = n(),
      M_subject = as.integer(all(M_row == 1L)),
      R_subject = as.integer(!any(cens_event_k == 1L)),
      p_denom = prod(p_obs),
      g_denom = prod(g_hat),
      Y_tau = min(max(Bk + Xk), tau),
      .groups = "drop"
    ) %>%
    mutate(
      W = ifelse(
        M_subject == 1L & R_subject == 1L,
        1 / pmax(p_denom * g_denom, 1e-8),
        0
      )
    )
}

# Estimates fixed-policy RMST up to tau using a normalized IPCW estimator.
estimate_ipcw_rmst_one_partition <- function(
    long_dat,
    policy_fun,
    history_vars = c("Z0", "Zk", "Ak_1", "Bk"),
    tau,
    num.trees.censor = 300,
    min.node.size = 10
) {
  validate_ipcw_data(long_dat, history_vars)

  dat <- long_dat %>% arrange(id, k)
  dat <- add_policy_match_ipcw(dat, policy_fun)

  prop_obj <- estimate_propensity_ipcw(dat, history_vars)
  dat$p_obs <- prop_obj$p_obs

  dat <- define_censoring_event_ipcw(dat, tau)
  censor_fit <- fit_censoring_forest_ipcw(
    dat,
    history_vars = history_vars,
    num.trees = num.trees.censor,
    min.node.size = min.node.size
  )
  dat$g_hat <- predict_censor_survival_ipcw(censor_fit, dat)

  subject_level <- build_subject_level_ipcw(dat, tau)
  weight_sum <- sum(subject_level$W)
  if (weight_sum <= 0) {
    stop("All IPCW subject weights are zero; check policy matching and censoring support.")
  }

  value <- sum(subject_level$W * subject_level$Y_tau) / weight_sum
  positive_weights <- subject_level$W[subject_level$W > 0]

  out <- list(
    value = value,
    subject_level = subject_level,
    row_level = dat,
    propensity_fit = prop_obj$fit,
    censor_fit = censor_fit,
    ess = effective_sample_size_ipcw(subject_level$W),
    match_rate = mean(subject_level$M_subject),
    weight_summary = summary(positive_weights),
    censor_survival_summary = summary(dat$g_hat),
    treatment_prob_summary = summary(dat$p_obs),
    tau = tau,
    history_vars = history_vars
  )
  class(out) <- "ipcw_rmst_one_partition"
  out
}

# Sanity checks
# source("SimData.R")
# pi_always_treat <- function(df) rep(1L, nrow(df))
# 
# fit_ipcw <- estimate_ipcw_rmst_one_partition(
#   long_dat = long_dat,
#   policy_fun = pi_always_treat,
#   history_vars = c("Z0", "Zk", "Ak_1", "Bk"),
#   tau = 730
# )
# 
# fit_ipcw$value
# fit_ipcw$ess
# fit_ipcw$match_rate
