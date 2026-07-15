# Example run for the RMST-first single-partition Algorithm 1 implementation.

source("WithinPartitionAlg_RMST.R")
source("SimData.R")

sim <- simu(n = 200, K = 5, tau = 730, seed = 1)
long_dat <- sim$long


# Placeholder evaluation policy.
# Replace this with the learned policy from your training procedure.
pi_always_treat <- function(df) {
  rep(1L, nrow(df))
}


# Fit the single-partition RMST-first recursive forest.
alg1_rmst_fit <- fit_algorithm1_rmst_one_partition(
  long_dat = long_dat,
  policy_fun = pi_always_treat,
  history_vars = c("Z0", "Zk", "Ak_1", "Bk"),
  tau = 730,
  n_grid = 100,
  num.trees.1stub = 300,
  num.trees.curve = 150,
  num.trees.rmst = 150,
  min.node.size = 10
)


# Value estimate using matched baseline starts.
baseline_value <- value_algorithm1_rmst(alg1_rmst_fit, baseline_only = TRUE)
print(baseline_value$value)


# Predicted RMSTs at baseline.
baseline_dat <- alg1_rmst_fit$data %>% filter(k == 1L)
baseline_rmst_pred <- predict(alg1_rmst_fit, baseline_dat)
print(head(baseline_rmst_pred))


# Optional: inspect propagated survival curves from the final internal curve model.
baseline_curve_pred <- predict_curves_algorithm1_rmst(alg1_rmst_fit, baseline_dat)
print(dim(baseline_curve_pred))


# Optional: inspect which augmentation path was used at each stub length.
print(lapply(alg1_rmst_fit$augmentation_info, function(x) table(x$augmentation_type)))
