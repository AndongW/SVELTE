# Within-Partition Recursive Forest Refitting (RMST-first)
#
# Single-partition implementation of Algorithm 1 for L = 1.
# This file leaves WithinPartitionAlg.R untouched and uses a simpler public
# refit target: RMST of recursively augmented survival curves.
#
# Required columns in long_dat:
# id, k, Ak, Xk, delta_k, gamma_k, plus all history_vars
#
# Required packages:
# install.packages(c("dplyr", "ranger", "survival"))

library(dplyr)
library(ranger)
library(survival)

# Validation and utilities----

# Validates required columns and binary indicators for Algorithm 1 input data.
validate_algorithm1_data <- function(long_dat, history_vars) {
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

# Evaluates a right-continuous step survival curve at arbitrary times.
step_eval_surv <- function(time_grid, surv_vec, t) {
  if (length(time_grid) == 0L) {
    return(rep(1, length(t)))
  }

  idx <- findInterval(t, time_grid)
  out <- ifelse(idx <= 0L, 1, surv_vec[pmax(idx, 1L)])
  out[t < 0] <- 1
  pmin(pmax(out, 0), 1)
}

# Enforces basic survival-curve constraints on a numeric vector.
monotone_surv <- function(s) {
  s <- pmin(pmax(s, 0), 1)
  cummin(s)
}

# Builds the working time grid used for recursive curve construction.
make_time_grid <- function(dat, tau = NULL, n_grid = 100) {
  max_t <- if (is.null(tau)) max(dat$Xk, na.rm = TRUE) else tau
  sort(unique(c(0, seq(0, max_t, length.out = n_grid))))
}

# Computes RMST from a survival curve on the working grid.
rmst_from_curve <- function(grid, surv_vec, tau = max(grid, na.rm = TRUE)) {
  keep <- grid <= tau
  grid_sub <- grid[keep]
  surv_sub <- surv_vec[keep]

  if (length(grid_sub) < 2L) {
    return(0)
  }

  sum(diff(grid_sub) * head(surv_sub, -1))
}



# Propensity and matching weights----


# Fits a binary propensity model and returns observed-treatment probabilities.
estimate_propensity <- function(dat, history_vars) {
  f <- as.formula(paste("Ak ~", paste(history_vars, collapse = " + ")))
  fit <- glm(f, data = dat, family = binomial())
  p1 <- as.numeric(predict(fit, type = "response"))
  p1 <- pmin(pmax(p1, 1e-4), 1 - 1e-4)
  ifelse(dat$Ak == 1L, p1, 1 - p1)
}

# Adds policy-match indicators and inverse propensity matching weights.
add_matching_weights <- function(dat, policy_fun, history_vars) {
  dat$pi_hat <- as.integer(policy_fun(dat))
  dat$p_obs <- estimate_propensity(dat, history_vars)
  dat$M <- as.integer(dat$Ak == dat$pi_hat)
  dat$w_star <- dat$M / dat$p_obs
  dat
}



# Stage-1 forest ----


# Fits the weighted 1-stub survival forest on matched starts.
fit_survival_forest_1stub <- function(dat, history_vars, num.trees = 500, min.node.size = 10) {
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

# Predicts 1-stub survival curves from the initial survival forest.
predict_surv_1stub <- function(fit, newdata, grid) {
  pr <- predict(fit, data = newdata)
  base_grid <- pr$unique.death.times
  S <- pr$survival
  if (is.null(dim(S))) {
    S <- matrix(S, nrow = nrow(newdata), byrow = TRUE)
  }

  out <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(grid))
  for (i in seq_len(nrow(newdata))) {
    out[i, ] <- step_eval_surv(base_grid, S[i, ], grid)
  }
  out[, 1] <- 1
  t(apply(out, 1, monotone_surv))
}



# Internal curve models----


# Fits the internal per-timepoint curve model used to propagate recursion.
fit_internal_curve_model <- function(dat, S_mat, grid, history_vars,
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

# Predicts survival curves from the internal per-timepoint curve model.
predict_internal_curve_model <- function(fit, newdata) {
  out <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(fit$grid))
  nd <- newdata[, fit$history_vars, drop = FALSE]

  for (j in seq_along(fit$grid)) {
    out[, j] <- predict(fit$fits[[j]], data = nd)$predictions
  }

  out[, 1] <- 1
  t(apply(out, 1, monotone_surv))
}

# Dispatches curve prediction across the stage-1 forest and internal curve models.
predict_curve_model <- function(fit_obj, newdata, grid = NULL) {
  if (inherits(fit_obj, "ranger")) {
    return(predict_surv_1stub(fit_obj, newdata, grid))
  }
  if (is.list(fit_obj) && identical(fit_obj$type, "curve_model")) {
    return(predict_internal_curve_model(fit_obj, newdata))
  }

  stop("Unknown curve model object.")
}



# RMST forests----


# Fits the regression forest used for public RMST prediction at each recursion stage.
fit_rmst_forest <- function(dat, y_rmst, history_vars,
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

# Predicts scalar RMST targets from a fitted RMST forest.
predict_rmst_forest <- function(fit, newdata, history_vars) {
  nd <- newdata[, history_vars, drop = FALSE]
  as.numeric(predict(fit, data = nd)$predictions)
}



# Stub eligibility and recursion helpers ----


# Identifies matched starts with at least m consecutive observed visits.
get_valid_starts <- function(dat, m) {
  dat %>%
    arrange(id, k) %>%
    group_by(id) %>%
    mutate(has_m_stub = vapply(k, function(k0) all((k0:(k0 + m - 1)) %in% k), logical(1))) %>%
    ungroup() %>%
    filter(has_m_stub, M == 1L)
}

# Retrieves arbitrary visit rows needed for recursive augmentation.
get_target_rows <- function(starts, dat, target_k) {
  target_df <- starts %>% transmute(id = id, k_target = target_k)
  target_df %>% left_join(dat, by = c("id" = "id", "k_target" = "k"))
}

# Builds a lookup table for valid matched stub starts.
make_valid_start_lookup <- function(dat, m) {
  starts <- get_valid_starts(dat, m)
  paste(starts$id, starts$k, sep = "::")
}

# Finds the earliest later matched visit that starts a valid (m-1)-stub.
find_concat_start <- function(start_row, valid_prev_lookup, m) {
  later_visits <- seq.int(start_row$k + 1L, start_row$k + m - 1L)
  later_keys <- paste(start_row$id, later_visits, sep = "::")
  eligible <- later_visits[later_keys %in% valid_prev_lookup]

  if (length(eligible) == 0L) {
    NA_integer_
  } else {
    min(eligible)
  }
}

# Computes the elapsed observed duration from visit start k to visit start k'.
compute_concat_boundary <- function(start_row, dat, k_prime) {
  rows <- dat[dat$id == start_row$id & dat$k >= start_row$k & dat$k < k_prime, , drop = FALSE]
  sum(rows$Xk)
}

# Computes the observed duration of the first m-1 visits of the target m-stub.
compute_extension_boundary <- function(start_row, dat, m) {
  rows <- dat[
    dat$id == start_row$id & dat$k >= start_row$k & dat$k <= (start_row$k + m - 2L),
    ,
    drop = FALSE
  ]
  sum(rows$Xk)
}

# Splices a continuation curve onto the current (m-1)-stub curve.
construct_augmented_curve_one <- function(S_current, S_cont, boundary, grid, eps = 1e-8) {
  out <- numeric(length(grid))
  S_boundary <- max(step_eval_surv(grid, S_current, boundary), eps)

  for (j in seq_along(grid)) {
    t <- grid[j]
    if (t < boundary) {
      out[j] <- step_eval_surv(grid, S_current, t)
    } else {
      out[j] <- S_boundary * step_eval_surv(grid, S_cont, t - boundary)
    }
  }

  out[1] <- 1
  monotone_surv(out)
}

# Vectorizes augmented-curve construction over all valid starts at stub length m.
construct_augmented_curves <- function(starts, dat, prev_curve_model, forest_1stub, grid, m) {
  valid_prev_lookup <- make_valid_start_lookup(dat, m - 1L)
  concat_k <- vapply(
    seq_len(nrow(starts)),
    function(r) find_concat_start(starts[r, , drop = FALSE], valid_prev_lookup, m),
    integer(1)
  )
  use_concat <- !is.na(concat_k)

  S_current <- predict_curve_model(prev_curve_model, starts, grid)
  S_cont <- matrix(NA_real_, nrow = nrow(starts), ncol = length(grid))
  boundaries <- numeric(nrow(starts))
  augmentation_type <- ifelse(use_concat, "concatenation", "extension")

  if (any(use_concat)) {
    concat_rows <- get_target_rows(starts[use_concat, , drop = FALSE], dat, concat_k[use_concat])
    S_cont[use_concat, ] <- predict_curve_model(prev_curve_model, concat_rows, grid)
    boundaries[use_concat] <- vapply(
      which(use_concat),
      function(r) compute_concat_boundary(starts[r, , drop = FALSE], dat, concat_k[r]),
      numeric(1)
    )
  }

  if (any(!use_concat)) {
    ext_k <- starts$k[!use_concat] + m - 1L
    ext_rows <- get_target_rows(starts[!use_concat, , drop = FALSE], dat, ext_k)
    S_cont[!use_concat, ] <- predict_surv_1stub(forest_1stub, ext_rows, grid)
    boundaries[!use_concat] <- vapply(
      which(!use_concat),
      function(r) compute_extension_boundary(starts[r, , drop = FALSE], dat, m),
      numeric(1)
    )
  }

  S_aug <- matrix(NA_real_, nrow = nrow(starts), ncol = length(grid))
  for (r in seq_len(nrow(starts))) {
    S_aug[r, ] <- construct_augmented_curve_one(
      S_current = S_current[r, ],
      S_cont = S_cont[r, ],
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
      augmentation_type = augmentation_type,
      concat_k = ifelse(use_concat, concat_k, NA_integer_),
      boundary = boundaries
    )
  )
}



# Main fit / predict / value methods----


# Fits the full RMST-first recursive forest algorithm in the single-partition case.
fit_algorithm1_rmst_one_partition <- function(
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
  validate_algorithm1_data(long_dat, history_vars)

  dat <- long_dat %>% arrange(id, k)
  dat <- add_matching_weights(dat, policy_fun, history_vars)
  grid <- make_time_grid(dat, tau = tau, n_grid = n_grid)
  tau_use <- max(grid, na.rm = TRUE)

  J1 <- dat %>% filter(M == 1L)
  if (nrow(J1) == 0L) {
    stop("No matched 1-stubs. Check policy positivity or policy_fun.")
  }

  forest_1stub <- fit_survival_forest_1stub(
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

  S1 <- predict_surv_1stub(forest_1stub, J1, grid)
  augmented_curves[["1"]] <- S1
  augmentation_info[["1"]] <- data.frame(
    id = J1$id,
    k = J1$k,
    m = 1L,
    augmentation_type = "initial",
    concat_k = NA_integer_,
    boundary = J1$Xk
  )
  rmst_targets[["1"]] <- apply(S1, 1, function(s) rmst_from_curve(grid, s, tau = tau_use))

  max_K <- max(dat$k, na.rm = TRUE)

  for (m in 2:max_K) {
    starts_m <- get_valid_starts(dat, m)
    if (nrow(starts_m) == 0L) {
      break
    }

    aug_m <- construct_augmented_curves(
      starts = starts_m,
      dat = dat,
      prev_curve_model = curve_models[[m - 1L]],
      forest_1stub = forest_1stub,
      grid = grid,
      m = m
    )
    S_m <- aug_m$S_aug

    y_m <- apply(S_m, 1, function(s) rmst_from_curve(grid, s, tau = tau_use))

    curve_models[[m]] <- fit_internal_curve_model(
      dat = starts_m,
      S_mat = S_m,
      grid = grid,
      history_vars = history_vars,
      num.trees = num.trees.curve,
      min.node.size = min.node.size
    )

    rmst_forests[[m]] <- fit_rmst_forest(
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
    history_vars = history_vars
  )
  class(out) <- "algorithm1_rmst_one_partition"
  out
}

# Returns final-stage RMST predictions from a fitted recursive RMST object.
predict.algorithm1_rmst_one_partition <- function(object, newdata, ...) {
  if (object$K_max <= 1L || is.null(object$final_rmst_forest)) {
    S_hat <- predict_surv_1stub(object$forest_1stub, newdata, object$grid)
    return(as.numeric(apply(S_hat, 1, function(s) rmst_from_curve(object$grid, s, tau = object$tau))))
  }

  predict_rmst_forest(object$final_rmst_forest, newdata, object$history_vars)
}

# Exposes propagated survival-curve predictions from the final internal curve model.
predict_curves_algorithm1_rmst <- function(object, newdata) {
  predict_curve_model(object$final_curve_model, newdata, object$grid)
}

# Averages predicted RMST values over baseline starts or all matched starts.
value_algorithm1_rmst <- function(fit, baseline_only = TRUE) {
  dat <- fit$data
  eval_dat <- if (baseline_only) {
    dat %>% filter(k == 1L, M == 1L)
  } else {
    dat %>% filter(M == 1L)
  }

  pred <- predict(fit, eval_dat)
  list(value = mean(pred), values = pred, eval_data = eval_dat)
}
