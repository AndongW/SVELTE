# Within-Partition Recursive Forest Refitting (RMST-first, 1-stub extension)
#
# Alternative single-partition implementation of Algorithm 1 for L = 1.
# This version keeps the same overall recursion but uses a simpler Step 2:
# every m-stub is formed by extending the predicted (m-1)-stub with a
# predicted 1-stub from the initial forest.

library(dplyr)
library(ranger)
library(survival)

# Validation and utilities----

validate_algorithm1_data_1stubext <- function(long_dat, history_vars) {
  required_cols <- unique(c("id", "k", "Ak", "Xk", "delta_k", "gamma_k", history_vars))
  missing_cols <- setdiff(required_cols, names(long_dat))
  if (length(missing_cols) > 0L) {
    stop("long_dat is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!all(long_dat$delta_k %in% c(0, 1, 0L, 1L, NA))) {
    stop("delta_k must be coded as 0/1.")
  }
  if (!all(long_dat$gamma_k %in% c(0, 1, 0L, 1L, NA))) {
    stop("gamma_k must be coded as 0/1.")
  }

  invisible(TRUE)
}

step_eval_surv_1stubext <- function(time_grid, surv_vec, t) {
  if (length(time_grid) == 0L) {
    return(rep(1, length(t)))
  }

  idx <- findInterval(t, time_grid)
  out <- ifelse(idx <= 0L, 1, surv_vec[pmax(idx, 1L)])
  out[t < 0] <- 1
  pmin(pmax(out, 0), 1)
}

monotone_surv_1stubext <- function(s) {
  s <- pmin(pmax(s, 0), 1)
  cummin(s)
}

make_time_grid_1stubext <- function(dat, tau = NULL, n_grid = 100) {
  max_t <- if (is.null(tau)) max(dat$Xk, na.rm = TRUE) else tau
  sort(unique(c(0, seq(0, max_t, length.out = n_grid))))
}

rmst_from_curve_1stubext <- function(grid, surv_vec, tau = max(grid, na.rm = TRUE)) {
  keep <- grid <= tau
  grid_sub <- grid[keep]
  surv_sub <- surv_vec[keep]

  if (length(grid_sub) < 2L) {
    return(0)
  }

  sum(diff(grid_sub) * head(surv_sub, -1))
}

# Propensity and matching weights----

estimate_propensity_1stubext <- function(dat, history_vars) {
  f <- as.formula(paste("Ak ~", paste(history_vars, collapse = " + ")))
  fit <- glm(f, data = dat, family = binomial())
  p1 <- as.numeric(predict(fit, type = "response"))
  p1 <- pmin(pmax(p1, 1e-4), 1 - 1e-4)
  ifelse(dat$Ak == 1L, p1, 1 - p1)
}

add_matching_weights_1stubext <- function(dat, policy_fun, history_vars) {
  dat$pi_hat <- as.integer(policy_fun(dat))
  dat$p_obs <- estimate_propensity_1stubext(dat, history_vars)
  dat$M <- as.integer(dat$Ak == dat$pi_hat)
  dat$w_star <- dat$M / dat$p_obs
  dat
}

# Stage-1 forest ----

fit_survival_forest_1stub_1stubext <- function(dat, history_vars, num.trees = 500, min.node.size = 10) {
  f <- as.formula(paste0("Surv(Xk, delta_k) ~ ", paste(history_vars, collapse = " + ")))
  ranger(
    formula = f,
    data = dat,
    case.weights = dat$w_star,
    num.trees = num.trees,
    min.node.size = min.node.size,
    respect.unordered.factors = "order",
    seed = 1
  )
}

predict_surv_1stub_1stubext <- function(fit, newdata, grid) {
  pr <- predict(fit, data = newdata)
  base_grid <- pr$unique.death.times
  S <- pr$survival
  if (is.null(dim(S))) {
    S <- matrix(S, nrow = nrow(newdata), byrow = TRUE)
  }

  out <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(grid))
  for (i in seq_len(nrow(newdata))) {
    out[i, ] <- step_eval_surv_1stubext(base_grid, S[i, ], grid)
  }
  out[, 1] <- 1
  t(apply(out, 1, monotone_surv_1stubext))
}

# Internal curve models----

fit_internal_curve_model_1stubext <- function(dat, S_mat, grid, history_vars,
                                               num.trees = 300, min.node.size = 10) {
  fits <- vector("list", length(grid))
  names(fits) <- paste0("t", seq_along(grid))
  train_df <- dat[, history_vars, drop = FALSE]

  for (j in seq_along(grid)) {
    df_j <- cbind(train_df, y = S_mat[, j])
    fits[[j]] <- ranger(
      y ~ .,
      data = df_j,
      case.weights = dat$w_star,
      num.trees = num.trees,
      min.node.size = min.node.size,
      seed = 1000 + j
    )
  }

  list(type = "curve_model", fits = fits, grid = grid, history_vars = history_vars)
}

predict_internal_curve_model_1stubext <- function(fit, newdata) {
  out <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(fit$grid))
  nd <- newdata[, fit$history_vars, drop = FALSE]

  for (j in seq_along(fit$grid)) {
    out[, j] <- predict(fit$fits[[j]], data = nd)$predictions
  }

  out[, 1] <- 1
  t(apply(out, 1, monotone_surv_1stubext))
}

predict_curve_model_1stubext <- function(fit_obj, newdata, grid = NULL) {
  if (inherits(fit_obj, "ranger")) {
    return(predict_surv_1stub_1stubext(fit_obj, newdata, grid))
  }
  if (is.list(fit_obj) && identical(fit_obj$type, "curve_model")) {
    return(predict_internal_curve_model_1stubext(fit_obj, newdata))
  }

  stop("Unknown curve model object.")
}

# RMST forests----

fit_rmst_forest_1stubext <- function(dat, y_rmst, history_vars,
                                     num.trees = 300, min.node.size = 10) {
  df <- cbind(dat[, history_vars, drop = FALSE], y_rmst = y_rmst)
  ranger(
    y_rmst ~ .,
    data = df,
    case.weights = dat$w_star,
    num.trees = num.trees,
    min.node.size = min.node.size,
    seed = 2026
  )
}

predict_rmst_forest_1stubext <- function(fit, newdata, history_vars) {
  nd <- newdata[, history_vars, drop = FALSE]
  as.numeric(predict(fit, data = nd)$predictions)
}

# Stub eligibility and recursion helpers ----

get_valid_starts_1stubext <- function(dat, m) {
  dat %>%
    arrange(id, k) %>%
    group_by(id) %>%
    mutate(has_m_stub = vapply(k, function(k0) all((k0:(k0 + m - 1)) %in% k), logical(1))) %>%
    ungroup() %>%
    filter(has_m_stub, M == 1L)
}

get_target_rows_1stubext <- function(starts, dat, target_k) {
  target_df <- starts %>% transmute(id = id, k_target = target_k)
  target_df %>% left_join(dat, by = c("id" = "id", "k_target" = "k"))
}

compute_extension_boundary_1stubext <- function(start_row, dat, m) {
  rows <- dat[
    dat$id == start_row$id & dat$k >= start_row$k & dat$k <= (start_row$k + m - 2L),
    ,
    drop = FALSE
  ]
  sum(rows$Xk)
}

construct_augmented_curve_one_1stubext <- function(S_current, S_ext, boundary, grid, eps = 1e-8) {
  out <- numeric(length(grid))
  S_boundary <- max(step_eval_surv_1stubext(grid, S_current, boundary), eps)

  for (j in seq_along(grid)) {
    t <- grid[j]
    if (t < boundary) {
      out[j] <- step_eval_surv_1stubext(grid, S_current, t)
    } else {
      out[j] <- S_boundary * step_eval_surv_1stubext(grid, S_ext, t - boundary)
    }
  }

  out[1] <- 1
  monotone_surv_1stubext(out)
}

construct_augmented_curves_1stubext <- function(starts, dat, prev_curve_model, forest_1stub, grid, m) {
  S_current <- predict_curve_model_1stubext(prev_curve_model, starts, grid)
  ext_k <- starts$k + m - 1L
  ext_rows <- get_target_rows_1stubext(starts, dat, ext_k)
  S_ext <- predict_surv_1stub_1stubext(forest_1stub, ext_rows, grid)
  boundaries <- vapply(
    seq_len(nrow(starts)),
    function(r) compute_extension_boundary_1stubext(starts[r, , drop = FALSE], dat, m),
    numeric(1)
  )

  S_aug <- matrix(NA_real_, nrow = nrow(starts), ncol = length(grid))
  for (r in seq_len(nrow(starts))) {
    S_aug[r, ] <- construct_augmented_curve_one_1stubext(
      S_current = S_current[r, ],
      S_ext = S_ext[r, ],
      boundary = boundaries[r],
      grid = grid
    )
  }

  colnames(S_aug) <- paste0("t", seq_along(grid))
  list(
    S_aug = S_aug,
    augmentation = data.frame(
      id = starts$id,
      k = starts$k,
      m = m,
      augmentation_type = "1stub_extension",
      extension_k = ext_k,
      boundary = boundaries
    )
  )
}

# Main fit / predict / value methods----

fit_algorithm1_rmst_one_partition_1stubext <- function(
    long_dat,
    policy_fun,
    history_vars = c("Z0", "Zk", "Ak_1", "Bk"),
    tau = NULL,
    n_grid = 100,
    num.trees.1stub = 500,
    num.trees.curve = 300,
    num.trees.rmst = 300,
    min.node.size = 10
) {
  validate_algorithm1_data_1stubext(long_dat, history_vars)

  dat <- long_dat %>% arrange(id, k)
  dat <- add_matching_weights_1stubext(dat, policy_fun, history_vars)
  grid <- make_time_grid_1stubext(dat, tau = tau, n_grid = n_grid)
  tau_use <- max(grid, na.rm = TRUE)

  J1 <- dat %>% filter(M == 1L)
  if (nrow(J1) == 0L) {
    stop("No matched 1-stubs. Check policy positivity or policy_fun.")
  }

  forest_1stub <- fit_survival_forest_1stub_1stubext(
    J1,
    history_vars = history_vars,
    num.trees = num.trees.1stub,
    min.node.size = min.node.size
  )

  curve_models <- list()
  curve_models[[1]] <- forest_1stub

  rmst_forests <- list()
  augmented_curves <- list()
  augmentation_info <- list()
  rmst_targets <- list()
  fit_data <- list(`1` = J1)

  S1 <- predict_surv_1stub_1stubext(forest_1stub, J1, grid)
  augmented_curves[["1"]] <- S1
  augmentation_info[["1"]] <- data.frame(
    id = J1$id,
    k = J1$k,
    m = 1L,
    augmentation_type = "initial",
    extension_k = J1$k,
    boundary = J1$Xk
  )
  rmst_targets[["1"]] <- apply(S1, 1, function(s) rmst_from_curve_1stubext(grid, s, tau = tau_use))

  max_K <- max(dat$k, na.rm = TRUE)

  for (m in 2:max_K) {
    starts_m <- get_valid_starts_1stubext(dat, m)
    if (nrow(starts_m) == 0L) {
      break
    }

    aug_m <- construct_augmented_curves_1stubext(
      starts = starts_m,
      dat = dat,
      prev_curve_model = curve_models[[m - 1L]],
      forest_1stub = forest_1stub,
      grid = grid,
      m = m
    )
    S_m <- aug_m$S_aug

    y_m <- apply(S_m, 1, function(s) rmst_from_curve_1stubext(grid, s, tau = tau_use))

    curve_models[[m]] <- fit_internal_curve_model_1stubext(
      dat = starts_m,
      S_mat = S_m,
      grid = grid,
      history_vars = history_vars,
      num.trees = num.trees.curve,
      min.node.size = min.node.size
    )

    rmst_forests[[m]] <- fit_rmst_forest_1stubext(
      dat = starts_m,
      y_rmst = y_m,
      history_vars = history_vars,
      num.trees = num.trees.rmst,
      min.node.size = min.node.size
    )

    augmented_curves[[as.character(m)]] <- S_m
    augmentation_info[[as.character(m)]] <- aug_m$augmentation
    rmst_targets[[as.character(m)]] <- y_m
    fit_data[[as.character(m)]] <- starts_m
  }

  K_max <- length(curve_models)
  final_rmst_forest <- rmst_forests[[K_max]]

  out <- list(
    forest_1stub = forest_1stub,
    curve_models = curve_models,
    rmst_forests = rmst_forests,
    final_rmst_forest = final_rmst_forest,
    final_curve_model = curve_models[[K_max]],
    K_max = K_max,
    grid = grid,
    tau = tau_use,
    data = dat,
    fit_data = fit_data,
    augmented_curves = augmented_curves,
    augmentation_info = augmentation_info,
    rmst_targets = rmst_targets,
    history_vars = history_vars,
    step2_method = "1stub_extension"
  )
  class(out) <- "algorithm1_rmst_one_partition_1stubext"
  out
}

predict.algorithm1_rmst_one_partition_1stubext <- function(object, newdata, ...) {
  if (object$K_max <= 1L || is.null(object$final_rmst_forest)) {
    S_hat <- predict_surv_1stub_1stubext(object$forest_1stub, newdata, object$grid)
    return(as.numeric(apply(S_hat, 1, function(s) rmst_from_curve_1stubext(object$grid, s, tau = object$tau))))
  }

  predict_rmst_forest_1stubext(object$final_rmst_forest, newdata, object$history_vars)
}

predict_curves_algorithm1_rmst_1stubext <- function(object, newdata) {
  predict_curve_model_1stubext(object$final_curve_model, newdata, object$grid)
}

value_algorithm1_rmst_1stubext <- function(fit, baseline_only = TRUE) {
  dat <- fit$data
  eval_dat <- if (baseline_only) {
    dat %>% filter(k == 1L, M == 1L)
  } else {
    dat %>% filter(M == 1L)
  }

  pred <- predict(fit, eval_dat)
  list(value = mean(pred), values = pred, eval_data = eval_dat)
}
