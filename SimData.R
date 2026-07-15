# Simulation Parameters:
# n: Number of subjects
# K: Maximum number of decision stages 
# tau: End of study censoring time 
# boundaries: Vector of partition boundaries in days
# T0: Baseline mean time-to-failure (in days) when linear predictor = 0
# U0: Baseline mean time-to-next-visit
# C0: Baseline mean time-to-censoring

# Subject level (baseline):
# id: Subject identifier
# l: Partition index {1,...,L}, where L = length(boundaries) + 1
# Z0: Baseline covariate
# time: Total observed follow-up time (final stage start + final stage duration)
# failed: 1 = failure occurred; 0 = no failure
# censored: 1 = censored (including end of study); 0 = not censored

# Stage-wise rows:
# --- Index variables ---
# id: subject index
# k: stage index
# l: partition
# --- History variables H_k ---
# Z0: Baseline covariate
# Zk: Time-varying covariate at stage k
# Ak_1: Treatment from previous stage
# Bk: Stage k start time (elapsed time since baseline)
# --- Treatment variable ---
# Ak: Treatment at stage k
# --- Time variables ---
# Tk: Time-to-failure (from stage start)
# Uk: Time-to-next-visit (advancement time)
# Ck: Time-to-censoring
# Xk: Observed stage duration (minimum of Tk, Uk, Ck)
# --- Indicator variables ---
# delta_k: 1 = not censored within stage (failure or advance); 0 = censored
# gamma_k: 1 = failure before next visit; 0 = observe next visit
# study_limit_censor: 1 = administratively censored because max visit limit K was reached

# Multi-stage survival DTR simulation 

simu <- function(
    n = 1000,
    K = 3,
    tau = 730,
    boundaries = c(365),
    seed = 1,
    T0 = 2000,      # baseline mean time-to-failure when eta=0, 
    U0 = 240,      # baseline mean time-to-next-visit when eta=0
    C0 = 2400       # baseline mean time-to-censoring when eta=0 (>730 because we dont want too many censored)
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

  # Assign Bk to partition intervals (-Inf, b1], (b1, b2], ..., (b_{L-1}, Inf).
  l_of_B <- function(Bk) findInterval(Bk, vec = boundaries, rightmost.closed = TRUE) + 1L
  
  # Baseline (subject level) ----
  id <- seq_len(n)
  Z0 <- rnorm(n, 0, 1)
  
  # initialize time-varying covariate at k=1
  Zk <- 0.5 * Z0 + rnorm(n, 0, 1)
  
  # stage start times and eligibility
  Bk <- rep(0, n)
  eligible <- rep(TRUE, n)
  
  # previous treatment Ak_1 (A_{k-1}); define Ak_1=0 at baseline
  Ak_1 <- rep(0L, n)
  
  # storage for long data
  out <- vector("list", n * K)
  idx <- 1L
  
  for (k in seq_len(K)) {
    at_risk <- which(eligible & (Bk < tau))
    if (length(at_risk) == 0) break
    
    i  <- at_risk
    z0 <- Z0[i]
    z  <- Zk[i]
    b  <- Bk[i]
    a_prev <- Ak_1[i]
    # time-based partition at stage start
    l_k <- l_of_B(b)
    
    # Treatment Ak ~ P(Ak=1 | Hk) ----
    # Hk = (Z0, Zk, Ak_1, Bk) in this simple simulation
    eta_A <- -0.2 + 0.5 * z0 + 0.6 * z + 0.2 * a_prev + 0.001 * b 
    # + 0.15 * (l_k - 1)
    pA    <- expit(eta_A)
    Ak    <- rbinom(length(i), 1, pA)
    
    # Latent times: Tk, Uk, Ck ----
    # We construct eta on log-rate scale; then scale = 1/exp(eta).
    # Failure hazard decreases with Ak and has heterogeneity via Ak:Zk.
    eta_T <- log(1 / T0) + 0.15 * z0 + 0.15 * z - 0.35 * Ak - 0.10 * Ak * z
    # + 0.10 * (l_k - 1)
    rate_T  <- exp(eta_T)
    Tk <- rexp(length(i), rate = rate_T)  
    
    # Next-visit: treated followed sooner (higher hazard => smaller mean)
    eta_U <- log(1 / U0) + 0.05 * z0 + 0.05 * z + 0.15 * Ak
    # + 0.05 * (l_k - 1)
    shape_U <- 3
    mu_U <- exp(-eta_U)
    scale_U <- mu_U / gamma(1 + 1 / shape_U)
    Uk <- rweibull(length(i), shape = shape_U, scale = scale_U)
    
    # Censoring: depends on Zk and partition (dependent censoring)
    eta_C <- log(1 / C0) + 0.05 * z0 + 0.10 * z + 0.00 * Ak
    # + 0.35 * (l_k - 1)
    rate_C  <- exp(eta_C)
    Ck <- rexp(length(i), rate = rate_C)
    
    # ---- Observed stage duration with end-of-study censoring ----
    rem_admin <- pmax(tau - b, 0)                 # remaining time to tau from stage start
    Xk <- pmin(Tk, Uk, Ck, rem_admin)         # observed stage duration (includes admin)
    
    # delta_k: 1 if NOT censored within stage (failure or advance),
    #          0 if censored (includes end-of-study admin censoring)
    admin_first <- (Xk >= rem_admin) & (rem_admin <= Tk) & (rem_admin <= Uk) & (rem_admin <= Ck)
    censor_first <- (!admin_first) & (Ck <= Tk) & (Ck <= Uk)
    delta_k <- as.integer(!admin_first & !censor_first)
    
    # gamma_k: 1 if failure before next visit, 0 otherwise (defined via latent Tk,Uk)
    gamma_k <- as.integer(Tk <= Uk)
    
    # stage advancement happens when Uk is the minimum and not censored/admin
    # equivalently: delta_k=1 and gamma_k=0 and Uk <= Tk (true) and Uk <= Ck and Uk < rem_admin
    advance <- (delta_k == 1L) & (gamma_k == 0L) & (Uk <= Ck) & (Uk < rem_admin)

    # If the subject would advance beyond the last allowed visit, treat this as
    # design-based censoring due to the study visit limit.
    study_limit_censor <- as.integer((k == K) & advance)
    if (any(study_limit_censor == 1L)) {
      delta_k[study_limit_censor == 1L] <- 0L
      gamma_k[study_limit_censor == 1L] <- 0L
      advance[study_limit_censor == 1L] <- FALSE
    }
    
    # update stage start time
    B_next <- b + Xk
    
    # update eligibility for next stage
    eligible[i] <- advance & (B_next < tau)
    
    # update Zk for those who advance
    Z_next <- 0.6 * z + 0.25 * z0 + 0.35 * Ak + rnorm(length(i), 0, 1)
    Zk[i] <- ifelse(eligible[i], Z_next, Zk[i])
    
    # update Ak_1 and Bk for those who advance
    Ak_1[i] <- ifelse(eligible[i], Ak, Ak_1[i])
    Bk[i] <- ifelse(eligible[i], B_next, Bk[i])
    
    # Write long rows ----
    for (j in seq_along(i)) {
      out[[idx]] <- data.frame(
        id       = i[j],
        k        = k,
        l        = l_k[j],
        # History Hk components
        Z0       = z0[j],
        Zk       = z[j],
        Ak_1     = a_prev[j],
        Bk       = b[j],
        # Treatment
        Ak       = Ak[j],
        # Time variables (latent + observed)
        Tk       = Tk[j],
        Uk       = Uk[j],
        Ck       = Ck[j],
        Xk       = Xk[j],
        # Indicators
        delta_k  = delta_k[j],
        gamma_k = gamma_k[j],
        study_limit_censor = study_limit_censor[j]
      )
      idx <- idx + 1L
    }
  }
  
  long_dat <- do.call(rbind, out[seq_len(idx - 1L)])
  long_dat <- long_dat[order(long_dat$id, long_dat$k), ]
  rownames(long_dat) <- NULL
  
  # Subject-level summary ----
  last_row <- long_dat[!duplicated(long_dat$id, fromLast = TRUE), ]
  
  time <- last_row$Bk + last_row$Xk
  
  # failure occurs if delta_k=1 and gamma_k=1 and Tk is the min among Tk,Uk,Ck,rem_admin
  rem_admin_last <- pmax(tau - last_row$Bk, 0)
  failed <- as.integer(
    (last_row$delta_k == 1L) &
      (last_row$gamma_k == 1L) &
      (last_row$Tk <= last_row$Ck) &
      (last_row$Tk <  rem_admin_last)
  )
  
  censored <- as.integer(failed == 0L)  # includes end-of-study censoring at tau and non-admin censoring
  
  subject_dat <- data.frame(
    id       = last_row$id,
    l        = last_row$l,
    Z0       = last_row$Z0,
    time     = time,
    failed   = failed,
    censored = censored
  )
  
  list(
    long = long_dat,
    subject = subject_dat,
    tau = tau,
    boundaries = boundaries,
    L = L
  )
}

plot_patient_trajectories <- function(sim, ids = 1:10, tau = NULL, boundaries = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use plot_patient_trajectories().")
  }

  long_dat <- sim$long
  if (is.null(tau)) {
    tau <- sim$tau
  }
  if (is.null(boundaries)) {
    boundaries <- sim$boundaries
  }

  if (is.null(long_dat)) {
    stop("sim must contain a long data frame in sim$long.")
  }

  if (is.null(tau)) {
    tau <- max(long_dat$Bk + long_dat$Xk)
    warning("tau was not found in sim; using max(Bk + Xk) from sim$long.")
  }

  if (is.null(boundaries)) {
    boundaries <- numeric(0)
    warning("boundaries were not found in sim; plotting without partition boundary lines.")
  }

  dat <- long_dat[long_dat$id %in% ids, ]
  dat <- dat[order(dat$id, dat$k), ]

  if (!"study_limit_censor" %in% names(dat)) {
    dat$study_limit_censor <- 0L
  }

  if (nrow(dat) == 0) {
    stop("No selected ids were found in sim$long.")
  }

  rem_admin <- pmax(tau - dat$Bk, 0)
  admin_first <- (dat$Xk >= rem_admin) & (rem_admin <= dat$Tk) &
    (rem_admin <= dat$Uk) & (rem_admin <= dat$Ck)
  censor_first <- (!admin_first) & (dat$Ck <= dat$Tk) & (dat$Ck <= dat$Uk)
  failure_first <- (dat$delta_k == 1L) & (dat$gamma_k == 1L) &
    (dat$Tk <= dat$Ck) & (dat$Tk < rem_admin)

  dat$x_start <- dat$Bk
  dat$x_end <- dat$Bk + dat$Xk
  dat$visit_label <- paste0("k=", dat$k)
  dat$end_type <- ifelse(
    failure_first,
    "Failure",
    ifelse(censor_first | admin_first | (dat$study_limit_censor == 1L), "Censored", "Ongoing")
  )
  dat$id_factor <- factor(dat$id, levels = rev(sort(unique(dat$id))))

  p <- ggplot2::ggplot(dat) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = x_start,
        xend = x_end,
        y = id_factor,
        yend = id_factor
      ),
      linewidth = 1.1,
      color = "grey35"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = x_start, y = id_factor, fill = factor(Ak)),
      shape = 21,
      size = 3,
      color = "black"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_start, y = id_factor, label = visit_label),
      nudge_y = 0.25,
      size = 3
    ) +
    # ggplot2::geom_point(
    #   data = dat[dat$end_type == "Advance", ],
    #   ggplot2::aes(x = x_end, y = id_factor),
    #   shape = 17,
    #   size = 2.8,
    #   color = "#1b9e77"
    # ) +
    ggplot2::geom_point(
      data = dat[dat$end_type == "Failure", ],
      ggplot2::aes(x = x_end, y = id_factor),
      shape = 4,
      stroke = 1.2,
      size = 3,
      color = "#d95f02"
    ) +
    ggplot2::geom_point(
      data = dat[dat$end_type == "Censored", ],
      ggplot2::aes(x = x_end, y = id_factor),
      shape = 1,
      stroke = 1.1,
      size = 3,
      color = "#7570b3"
    ) +
    # ggplot2::geom_vline(
    #   xintercept = boundaries,
    #   linetype = "dashed",
    #   color = "grey55"
    # ) +
    ggplot2::geom_vline(
      xintercept = tau,
      linetype = "dotted",
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = c("0" = "white", "1" = "black"),
      name = "Treatment Ak",
      labels = c("0" = "No", "1" = "Yes")
    ) +
    ggplot2::labs(
      x = "Time",
      y = "Patient",
      title = "Patient Trajectories"
    ) +
    ggplot2::theme_bw()

  if (length(boundaries) > 0) {
    partition_df <- data.frame(
      x = c(0, boundaries),
      xend = c(boundaries, tau),
      label = paste("Partition", seq_along(c(0, boundaries)))
    )
    partition_df$xmid <- (partition_df$x + partition_df$xend) / 2

    p <- p + ggplot2::annotate(
      "text",
      x = partition_df$xmid,
      y = length(unique(dat$id)) + 0.8,
      label = partition_df$label,
      size = 3
    )
  }

  p
}

summarize_sim_diagnostics <- function(sim, tau = NULL) {
  long_dat <- sim$long
  subject_dat <- sim$subject

  if (is.null(long_dat) || is.null(subject_dat)) {
    stop("sim must contain sim$long and sim$subject.")
  }

  if (is.null(tau)) {
    tau <- sim$tau
  }

  if (is.null(tau)) {
    tau <- max(long_dat$Bk + long_dat$Xk)
    warning("tau was not found in sim; using max(Bk + Xk) from sim$long.")
  }

  n_subjects <- nrow(subject_dat)
  visit_counts <- as.integer(table(long_dat$id))
  visit_count_props <- prop.table(table(visit_counts))

  last_rows <- long_dat[!duplicated(long_dat$id, fromLast = TRUE), ]
  if (!"study_limit_censor" %in% names(last_rows)) {
    last_rows$study_limit_censor <- 0L
  }
  early_terminal <- mean(last_rows$k == 1L)
  reached_tau <- mean(abs(subject_dat$time - tau) < 1e-8)
  study_limit_censor <- mean(last_rows$study_limit_censor == 1L)

  list(
    n_subjects = n_subjects,
    visit_count_table = table(visit_counts),
    visit_count_proportions = visit_count_props,
    proportion_two_or_three_visits = mean(visit_counts %in% c(2L, 3L)),
    proportion_one_visit = mean(visit_counts == 1L),
    proportion_terminal_on_visit_1 = early_terminal,
    proportion_reaching_tau = reached_tau,
    proportion_study_limit_censored = study_limit_censor
  )
}

# ---- Example ----
# n = 1000
# K = 3
# tau = 730
# boundaries = c(365)
# seed = 1990
# T0 = 600
# U0 = 180
# C0 = 900
sim1 <- simu(boundaries = c(365), seed = 1990)

long_dat    <- sim1$long
subject_dat <- sim1$subject

# sanity checks
# table(long_dat$k)
# with(subject_dat, table(failed, censored))
# head(long_dat, 10)
# head(subject_dat, 10)
# hist(subject_dat$time)

# Diagnostics for monitoring the simulated population profile
sim_diag <- summarize_sim_diagnostics(sim1)
sim_diag$visit_count_table
round(sim_diag$visit_count_proportions, 3)
sim_diag$proportion_two_or_three_visits
sim_diag$proportion_one_visit
sim_diag$proportion_terminal_on_visit_1
sim_diag$proportion_reaching_tau
sim_diag$proportion_study_limit_censored

# Example plot for the first 10 patients
plot_patient_trajectories(sim1, ids = 1:20)
