# Comparison workflow utilities for SVELTE RMST recursion versus IPCW RMST.

source("WithinPartitionAlg_RMST_1StubExt.R")
source("IPCW_RMST.R")

library(dplyr)


# Loads the simulation generator only when replicate workflows need it.
ensure_simu_available_svelte <- function() {
  if (!exists("simu", mode = "function")) {
    source("SimData.R")
  }
  invisible(TRUE)
}


# Infers the default history variables from the available baseline column names.
infer_history_vars_svelte <- function(long_dat) {
  if (all(c("X0", "Zk", "Ak_1", "Bk") %in% names(long_dat))) {
    return(c("X0", "Zk", "Ak_1", "Bk"))
  }
  if (all(c("Z0", "Zk", "Ak_1", "Bk") %in% names(long_dat))) {
    return(c("Z0", "Zk", "Ak_1", "Bk"))
  }
  stop("Could not infer history_vars from long_dat.")
}

# Simulates long-format trajectories under a fixed policy using the current data-generating mechanism.
simulate_fixed_policy_data <- function(
    n = 1000,
    K = 3,
    tau = 730,
    boundaries = c(365),
    seed = 1,
    T0 = 2000,
    U0 = 240,
    C0 = 2400,
    policy_fun,
    censoring = TRUE
) {
  set.seed(seed)
  expit <- function(x) 1 / (1 + exp(-x))

  if (!is.numeric(boundaries) || any(!is.finite(boundaries))) {
    stop("boundaries must be a numeric vector of finite values.")
  }

  boundaries <- sort(unique(boundaries))
  if (any(boundaries <= 0) || any(boundaries >= tau)) {
    stop("All boundaries must be strictly between 0 and tau.")
  }

  L <- length(boundaries) + 1L
  l_of_B <- function(Bk) findInterval(Bk, vec = boundaries, rightmost.closed = TRUE) + 1L

  id <- seq_len(n)
  Z0 <- rnorm(n, 0, 1)
  Zk <- 0.5 * Z0 + rnorm(n, 0, 1)
  Bk <- rep(0, n)
  eligible <- rep(TRUE, n)
  Ak_1 <- rep(0L, n)

  out <- vector("list", n * K)
  idx <- 1L

  for (k in seq_len(K)) {
    at_risk <- which(eligible & (Bk < tau))
    if (length(at_risk) == 0L) {
      break
    }

    i <- at_risk
    z0 <- Z0[i]
    z <- Zk[i]
    b <- Bk[i]
    a_prev <- Ak_1[i]
    l_k <- l_of_B(b)

    hist_df <- data.frame(
      id = i,
      k = k,
      l = l_k,
      Z0 = z0,
      Zk = z,
      Ak_1 = a_prev,
      Bk = b
    )
    Ak <- as.integer(policy_fun(hist_df))

    eta_T <- log(1 / T0) + 0.15 * z0 + 0.15 * z - 0.35 * Ak - 0.10 * Ak * z
    rate_T <- exp(eta_T)
    Tk <- rexp(length(i), rate = rate_T)

    eta_U <- log(1 / U0) + 0.05 * z0 + 0.05 * z + 0.15 * Ak
    shape_U <- 3
    mu_U <- exp(-eta_U)
    scale_U <- mu_U / gamma(1 + 1 / shape_U)
    Uk <- rweibull(length(i), shape = shape_U, scale = scale_U)

    if (censoring) {
      eta_C <- log(1 / C0) + 0.05 * z0 + 0.10 * z + 0.00 * Ak
      rate_C <- exp(eta_C)
      Ck <- rexp(length(i), rate = rate_C)
    } else {
      Ck <- rep(Inf, length(i))
    }

    rem_admin <- pmax(tau - b, 0)
    Xk <- pmin(Tk, Uk, Ck, rem_admin)

    admin_first <- (Xk >= rem_admin) & (rem_admin <= Tk) & (rem_admin <= Uk) & (rem_admin <= Ck)
    censor_first <- (!admin_first) & (Ck <= Tk) & (Ck <= Uk)
    delta_k <- as.integer(!admin_first & !censor_first)
    gamma_k <- as.integer(Tk <= Uk)
    advance <- (delta_k == 1L) & (gamma_k == 0L) & (Uk <= Ck) & (Uk < rem_admin)

    study_limit_censor <- as.integer((k == K) & advance)
    if (any(study_limit_censor == 1L)) {
      delta_k[study_limit_censor == 1L] <- 0L
      gamma_k[study_limit_censor == 1L] <- 0L
      advance[study_limit_censor == 1L] <- FALSE
    }

    B_next <- b + Xk
    eligible[i] <- advance & (B_next < tau)

    Z_next <- 0.6 * z + 0.25 * z0 + 0.35 * Ak + rnorm(length(i), 0, 1)
    Zk[i] <- ifelse(eligible[i], Z_next, Zk[i])
    Ak_1[i] <- ifelse(eligible[i], Ak, Ak_1[i])
    Bk[i] <- ifelse(eligible[i], B_next, Bk[i])

    for (j in seq_along(i)) {
      out[[idx]] <- data.frame(
        id = i[j],
        k = k,
        l = l_k[j],
        Z0 = z0[j],
        Zk = z[j],
        Ak_1 = a_prev[j],
        Bk = b[j],
        Ak = Ak[j],
        Tk = Tk[j],
        Uk = Uk[j],
        Ck = Ck[j],
        Xk = Xk[j],
        delta_k = delta_k[j],
        gamma_k = gamma_k[j],
        study_limit_censor = study_limit_censor[j]
      )
      idx <- idx + 1L
    }
  }

  long_dat <- do.call(rbind, out[seq_len(idx - 1L)])
  long_dat <- long_dat[order(long_dat$id, long_dat$k), ]
  rownames(long_dat) <- NULL

  last_row <- long_dat[!duplicated(long_dat$id, fromLast = TRUE), ]
  time <- pmin(last_row$Bk + last_row$Xk, tau)
  rem_admin_last <- pmax(tau - last_row$Bk, 0)
  failed <- as.integer(
    (last_row$delta_k == 1L) &
      (last_row$gamma_k == 1L) &
      (last_row$Tk <= last_row$Ck) &
      (last_row$Tk < rem_admin_last)
  )
  censored_subject <- as.integer(failed == 0L)

  subject_dat <- data.frame(
    id = last_row$id,
    l = last_row$l,
    Z0 = last_row$Z0,
    time = time,
    failed = failed,
    censored = censored_subject
  )

  list(
    long = long_dat,
    subject = subject_dat,
    tau = tau,
    boundaries = boundaries,
    L = L
  )
}

# Fits the recursive RMST estimator and averages predictions over all baseline subjects.
estimate_svelte_rmst_value <- function(
    long_dat,
    policy_fun,
    history_vars,
    tau,
    n_grid = 100,
    num.trees.1stub = 300,
    num.trees.curve = 150,
    num.trees.rmst = 150,
    min.node.size = 10
) {
  fit <- fit_algorithm1_rmst_one_partition_1stubext(
    long_dat = long_dat,
    policy_fun = policy_fun,
    history_vars = history_vars,
    tau = tau,
    n_grid = n_grid,
    num.trees.1stub = num.trees.1stub,
    num.trees.curve = num.trees.curve,
    num.trees.rmst = num.trees.rmst,
    min.node.size = min.node.size
  )

  baseline_dat <- fit$data %>% filter(k == 1L)
  baseline_pred <- predict(fit, baseline_dat)

  list(
    fit = fit,
    value = mean(baseline_pred),
    baseline_predictions = baseline_pred,
    baseline_data = baseline_dat
  )
}

# Computes the uncensored truth for a fixed policy by large-sample policy simulation.
estimate_true_policy_rmst <- function(
    policy_fun,
    n_truth = 10000,
    K = 3,
    tau = 730,
    boundaries = c(365),
    seed = 999,
    T0 = 2000,
    U0 = 240
) {
  sim_truth <- simulate_fixed_policy_data(
    n = n_truth,
    K = K,
    tau = tau,
    boundaries = boundaries,
    seed = seed,
    T0 = T0,
    U0 = U0,
    C0 = Inf,
    policy_fun = policy_fun,
    censoring = FALSE
  )

  list(
    value = mean(pmin(sim_truth$subject$time, tau)),
    sim = sim_truth
  )
}

# Runs one replicate of the SVELTE-versus-IPCW RMST comparison for a fixed policy.
run_svelte_ipcw_comparison_once <- function(
    long_dat,
    policy_fun,
    tau,
    history_vars = NULL,
    truth_args = NULL,
    svelte_args = list(),
    ipcw_args = list()
) {
  if (is.null(history_vars)) {
    history_vars <- infer_history_vars_svelte(long_dat)
  }

  svelte_defaults <- list(
    long_dat = long_dat,
    policy_fun = policy_fun,
    history_vars = history_vars,
    tau = tau
  )
  svelte_fit <- do.call(estimate_svelte_rmst_value, modifyList(svelte_defaults, svelte_args))

  ipcw_defaults <- list(
    long_dat = long_dat,
    policy_fun = policy_fun,
    history_vars = history_vars,
    tau = tau
  )
  ipcw_fit <- do.call(estimate_ipcw_rmst_one_partition, modifyList(ipcw_defaults, ipcw_args))

  truth <- NULL
  if (!is.null(truth_args)) {
    truth_defaults <- list(policy_fun = policy_fun, tau = tau)
    truth <- do.call(estimate_true_policy_rmst, modifyList(truth_defaults, truth_args))
  }

  out <- list(
    svelte = svelte_fit,
    ipcw = ipcw_fit,
    truth = truth,
    summary = data.frame(
      svelte_value = svelte_fit$value,
      ipcw_value = ipcw_fit$value,
      truth_value = if (is.null(truth)) NA_real_ else truth$value,
      ipcw_ess = ipcw_fit$ess,
      ipcw_match_rate = ipcw_fit$match_rate,
      svelte_K_max = svelte_fit$fit$K_max
    )
  )
  class(out) <- "svelte_ipcw_comparison"
  out
}

# Repeats the fixed-policy comparison across multiple observed-data replicates.
run_svelte_ipcw_comparison_replicates <- function(
    n_reps,
    policy_fun,
    sim_args,
    history_vars = NULL,
    truth_args = NULL,
    svelte_args = list(),
    ipcw_args = list(),
    seed_offset = 0
) {
  ensure_simu_available_svelte()

  results <- vector("list", n_reps)

  for (r in seq_len(n_reps)) {
    sim_r <- do.call(simu, modifyList(sim_args, list(seed = seed_offset + r)))
    results[[r]] <- run_svelte_ipcw_comparison_once(
      long_dat = sim_r$long,
      policy_fun = policy_fun,
      tau = sim_r$tau,
      history_vars = history_vars,
      truth_args = truth_args,
      svelte_args = svelte_args,
      ipcw_args = ipcw_args
    )
  }

  summaries <- bind_rows(lapply(results, function(x) x$summary), .id = "rep")
  list(results = results, summary = summaries)
}


# Example fixed-policy comparison run.
ensure_simu_available_svelte()

pi_always_treat <- function(df) {
  rep(1L, nrow(df))
}

sim_comp_example <- simu(
  n = 200,
  K = 3,
  tau = 730,
  boundaries = c(365),
  seed = 2026
)
# 
# comp_example <- run_svelte_ipcw_comparison_once(
#   long_dat = sim_comp_example$long,
#   policy_fun = pi_always_treat,
#   tau = sim_comp_example$tau,
#   history_vars = c("Z0", "Zk", "Ak_1", "Bk"),
#   truth_args = list(
#     n_truth = 2000,
#     K = 3,
#     tau = 730,
#     boundaries = c(365),
#     seed = 4040
#   ),
#   svelte_args = list(
#     n_grid = 60,
#     num.trees.1stub = 80,
#     num.trees.curve = 40,
#     num.trees.rmst = 40,
#     min.node.size = 5
#   ),
#   ipcw_args = list(
#     num.trees.censor = 80,
#     min.node.size = 5
#   )
# )
# 
# print(comp_example$summary)

res <- run_svelte_ipcw_comparison_replicates(
  n_reps = 50,
  policy_fun = pi_always_treat,
  sim_args = list(
    n = 200,
    K = 3,
    tau = 730,
    boundaries = c(365)
  ),
  history_vars = c("Z0", "Zk", "Ak_1", "Bk"),
  truth_args = list(
    n_truth = 2000,
    K = 3,
    tau = 730,
    boundaries = c(365),
    seed = 4040
  ),
  svelte_args = list(
    n_grid = 60,
    num.trees.1stub = 80,
    num.trees.curve = 40,
    num.trees.rmst = 40,
    min.node.size = 5
  ),
  ipcw_args = list(
    num.trees.censor = 80,
    min.node.size = 5
  ),
  seed_offset = 5000
)

within(res$summary, {
  svelte_err <- svelte_value - truth_value
  ipcw_err <- ipcw_value - truth_value
})

transform(
  data.frame(
    method = c("SVELTE", "IPCW"),
    bias = c(mean(res$summary$svelte_value - res$summary$truth_value),
             mean(res$summary$ipcw_value - res$summary$truth_value)),
    rmse = c(sqrt(mean((res$summary$svelte_value - res$summary$truth_value)^2)),
             sqrt(mean((res$summary$ipcw_value - res$summary$truth_value)^2))),
    sd = c(sd(res$summary$svelte_value),
           sd(res$summary$ipcw_value))
  )
)
