# LogisticSDEtools-R-
R code to simulate a logistic stochastic differential equation using Euler–Maruyama in log-space. Provides parameter estimation routines, Monte Carlo experiments, and visualization tools. Fully self-contained, using only base R.

# logisticSDEtools

R code to simulate a logistic stochastic differential equation using Euler–Maruyama in log-space.  
Provides parameter estimation routines, Monte Carlo experiments, and visualization tools.  
Fully self-contained, using only base R.

## Features
- Simulation of logistic SDE trajectories.
- Parameter estimation from simulated data.
- Monte Carlo experiments to study estimator behavior.
- Visualization of results (base R plots).

## Quick start
# =========================================================
# Logistic SDE — Simulation (Euler–Maruyama) and Point Estimation
#
# MODEL (46):
#   dY_t = [r - eta * Y_t] * Y_t dt + sigma * Y_t dB_t,   Y(0) > 0.
#
# INITIAL CONDITION (47): Stationary Gamma
#   Y_0 ~ Gamma(u, v)  with  u = 2r/sigma^2 - 1  (shape),  v = 2*eta/sigma^2  (rate).
#   Note: we use the "rate" parameterization (mean = u / v).
#
# ESTIMATOR OF sigma (50) using Z_t = log(Y_t):
#   hat_sigma^2 = (1 / (N * dt)) * sum_{n=0}^{N-1} ( Z_{n+1} - Z_n )^2
#
# ESTIMATOR OF eta (53):
#   hat_eta = [ ( log(Y_T / Y_0)/T ) * ∫_0^T Y(s) ds  -  ∫_0^T Y(s) dZ(s) ]
#             / [ ∫_0^T Y(s)^2 ds  - (1/T) * ( ∫_0^T Y(s) ds )^2 ]
#
# ESTIMATOR OF r (54):
#   hat_r = log(Y_T / Y_0)/T  +  hat_eta * ∫_0^T Y(s) ds  +  (sigma^2)/2
#
# NUMERICAL DISCRETIZATION used for (53)–(54):
#   ∫_0^T Y(s) ds     ≈  sum_{n=0}^{N-1} Y_{t_n} * dt
#   ∫_0^T Y(s)^2 ds   ≈  sum_{n=0}^{N-1} Y_{t_n}^2 * dt
#   ∫_0^T Y(s) dZ(s)  ≈  sum_{n=0}^{N-1} Y_{t_n} * (Z_{t_{n+1}} - Z_{t_n})   (Itô, pre-point)
#
# SIMULATION SCHEME (Euler–Maruyama in log-space to enforce Y_t > 0):
#   Define Z_t = log Y_t. Then
#     Z_{n+1} = Z_n + [ r - eta * exp(Z_n) - (sigma^2)/2 ] * dt + sigma * sqrt(dt) * xi_n,
#     Y_{n+1} = exp(Z_{n+1}),   with xi_n ~ N(0,1) i.i.d.
# =========================================================


# ---- Base packages check (they are part of base R; this is just a friendly guard) ----
if (!requireNamespace("stats", quietly = TRUE)) {
  stop("'stats' package should be available in base R.")
}
if (!requireNamespace("graphics", quietly = TRUE)) {
  stop("'graphics' package should be available in base R.")
}


# -----------------------------------------------------------------------------
# simulate_logistic_sde_em()
# Purpose: Simulate the logistic SDE (46) using Euler–Maruyama in log-space.
# Inputs:
#   r, eta, sigma : positive parameters of the model (46)
#   T, dt         : time horizon and time step; N := round(T/dt)
#   n_paths       : number of independent trajectories to simulate
#   seed          : optional RNG seed (integer) for reproducibility
#   y0            : optional initial value(s) > 0. If NULL and
#                   use_stationary_init=TRUE, draw from (47) Gamma(u,v).
#   use_stationary_init : if TRUE and y0 is NULL, sample Y0 ~ Gamma(u,v) per (47)
# Output:
#   A list with time grid, simulated Y and Z=log(Y), and metadata (params).
# Notes:
#   - Simulation uses Z to guarantee Y > 0.
#   - If using the stationary init (47), condition 2r/sigma^2 > 1 must hold.
# -----------------------------------------------------------------------------
simulate_logistic_sde_em <- function(r, eta, sigma, T, dt,
                                     n_paths = 1L,
                                     seed = NULL,
                                     y0 = NULL,
                                     use_stationary_init = TRUE) {
  stopifnot(r > 0, eta > 0, sigma > 0, T > 0, dt > 0, n_paths >= 1)
  if (!is.null(seed)) set.seed(seed)
  
  # Time grid
  N <- as.integer(round(T / dt))    # number of steps
  time <- seq(0, N) * dt            # length N+1 including t=0
  
  # Initial condition (47): Gamma(u, v) with 'rate' v
  Y0 <- numeric(n_paths)
  if (!is.null(y0)) {
    stopifnot(length(y0) %in% c(1, n_paths))
    Y0 <- if (length(y0) == 1) rep(y0, n_paths) else y0
    if (any(Y0 <= 0)) stop("All initial values y0 must be > 0.")
  } else if (use_stationary_init) {
    u <- 2 * r / (sigma^2) - 1   # shape
    v <- 2 * eta / (sigma^2)     # rate
    if (u <= 0) stop("Stationary init requires 2r/sigma^2 > 1 (i.e., u > 0).")
    Y0 <- stats::rgamma(n_paths, shape = u, rate = v)
  } else {
    Y0 <- rep(1.0, n_paths)      # neutral fallback
  }
  
  # Containers
  Y <- matrix(NA_real_, nrow = N + 1L, ncol = n_paths)  # levels
  Z <- matrix(NA_real_, nrow = N + 1L, ncol = n_paths)  # logs
  Y[1, ] <- Y0
  Z[1, ] <- log(Y0)
  
  # Euler–Maruyama in log-space (scheme detailed above)
  drift_const <- r - 0.5 * sigma^2
  sqdt <- sqrt(dt)
  
  for (n in 1:N) {
    xi <- stats::rnorm(n_paths)
    Z[n + 1, ] <- Z[n, ] + (drift_const - eta * exp(Z[n, ])) * dt + sigma * sqdt * xi
    Y[n + 1, ] <- exp(Z[n + 1, ])
  }
  
  list(
    time = time,
    Y = Y,
    Z = Z,
    params = list(r = r, eta = eta, sigma = sigma, T = N * dt, dt = dt, N = N, n_paths = n_paths)
  )
}


# -----------------------------------------------------------------------------
# estimate_sigma()
# Purpose: Estimate sigma via (50) using Z = log(Y).
# Equation (50):
#   hat_sigma^2 = (1 / (N * dt)) * sum_{n=0}^{N-1} ( Z_{n+1} - Z_n )^2
# Output:
#   Vector of length n_paths: hat_sigma for each path.
# -----------------------------------------------------------------------------
estimate_sigma <- function(sim) {
  Z  <- sim$Z
  dt <- sim$params$dt
  N  <- sim$params$N
  
  dZ <- Z[-1, , drop = FALSE] - Z[-nrow(Z), , drop = FALSE]
  sigma2_hat <- colSums(dZ^2) / (N * dt)
  sqrt(sigma2_hat)
}


# -----------------------------------------------------------------------------
# estimate_eta_eq53()
# Purpose: Estimate eta using (53).
# Discretization rules (see header) for ∫Y ds, ∫Y^2 ds and ∫Y dZ (Itô pre-point).
# Equation (53):
#   hat_eta = [ ( (Z_T - Z_0)/T ) * ∫_0^T Y(s) ds  -  ∫_0^T Y(s) dZ(s) ]
#             / [ ∫_0^T Y(s)^2 ds  - (1/T) * ( ∫_0^T Y(s) ds )^2 ].
# Output:
#   Vector of length n_paths: hat_eta for each path.
# -----------------------------------------------------------------------------
estimate_eta_eq53 <- function(sim) {
  Y  <- sim$Y
  Z  <- sim$Z
  dt <- sim$params$dt
  N  <- sim$params$N
  Tt <- N * dt
  
  Y_left <- Y[1:N, , drop = FALSE]                           # Y(t_n), n=0..N-1
  dZ     <- Z[2:(N+1), , drop = FALSE] - Z[1:N, , drop = FALSE]
  intY   <- colSums(Y_left)  * dt                            # ∫ Y ds
  intY2  <- colSums(Y_left^2) * dt                           # ∫ Y^2 ds
  intYdZ <- colSums(Y_left * dZ)                             # ∫ Y dZ
  Zdiff  <- Z[N+1, ] - Z[1, ]                                # Z_T - Z_0 = log(Y_T/Y_0)
  
  denom <- intY2 - (intY^2) / Tt
  if (any(denom <= 0)) warning("Non-positive denominator in (53) for some path(s).")
  
  num <- (Zdiff / Tt) * intY - intYdZ
  as.numeric(num / denom)
}


# -----------------------------------------------------------------------------
# estimate_r_eq54()
# Purpose: Estimate r using (54), plugging in hat_sigma from (50) and hat_eta from (53).
# Equation (54):
#   hat_r = (Z_T - Z_0)/T + hat_eta * ∫_0^T Y(s) ds + (hat_sigma^2)/2
# Output:
#   Vector of length n_paths: hat_r for each path.
# -----------------------------------------------------------------------------
estimate_r_eq54 <- function(sim, sigma_hat, eta_hat) {
  Z  <- sim$Z
  Y  <- sim$Y
  dt <- sim$params$dt
  N  <- sim$params$N
  Tt <- N * dt
  
  if (length(sigma_hat) == 1L) sigma_hat <- rep(sigma_hat, ncol(Y))
  if (length(eta_hat)   == 1L) eta_hat   <- rep(eta_hat,   ncol(Y))
  
  Y_left <- Y[1:N, , drop = FALSE]
  intY   <- colSums(Y_left) * dt
  Zdiff  <- Z[N+1, ] - Z[1, ]
  
  r_hat <- 0.5 * (sigma_hat^2) + (Zdiff + eta_hat * intY) / Tt
  as.numeric(r_hat)
}


# -----------------------------------------------------------------------------
# run_replication()
# Purpose: One end-to-end replication:
#   1) Simulate a single trajectory using (46) with init (47).
#   2) Estimate sigma via (50).
#   3) Estimate eta via (53).
#   4) Estimate r via (54) using the previous estimates.
# Output:
#   A named list with sigma_hat, eta_hat, r_hat.
# -----------------------------------------------------------------------------
run_replication <- function(r, eta, sigma, T, dt, seed = NULL) {
  sim <- simulate_logistic_sde_em(r = r, eta = eta, sigma = sigma,
                                  T = T, dt = dt, n_paths = 1L,
                                  seed = seed, use_stationary_init = TRUE)
  sigma_hat <- estimate_sigma(sim)[1]                   # (50)
  eta_hat   <- estimate_eta_eq53(sim)[1]                # (53)
  r_hat     <- estimate_r_eq54(sim, sigma_hat, eta_hat)[1]  # (54)
  
  list(sigma_hat = sigma_hat, eta_hat = eta_hat, r_hat = r_hat)
}


# --------------------------
# Minimal usage example
# --------------------------
# This small run simulates one path and returns the three point estimates.
# (Adjust T and dt as needed; larger T typically reduces variance.)
res1 <- run_replication(r = 2.0, eta = 0.5, sigma = 0.6,
                        T = 10, dt = 0.01, seed = 1)
res1


# =========================================================
# Step 2 — Monte Carlo with K replications (equation-cited)
#
# Goal:
#   Quantify the sampling variability of the point estimators derived in Step 1.
#   Each replication:
#     (i)  simulates one trajectory of the SDE (46) with initial law (47),
#     (ii) estimates sigma via (50),
#     (iii) estimates eta via (53),
#     (iv) estimates r via (54).
#   We then aggregate the K estimates to report mean, standard deviation,
#   and a simple “95% band” as mean ± 2·sd (dispersion of the estimator).
#
# References:
#   (46) SDE dynamics for Y_t
#   (47) Stationary Gamma initial condition for Y_0
#   (50) Estimator of sigma using Z = log Y
#   (53) Estimator of eta using discrete approximations of ∫Y ds, ∫Y^2 ds, ∫Y dZ
#   (54) Estimator of r plugging in hat_sigma and hat_eta
# =========================================================

# -----------------------------------------------------------------------------
# run_montecarlo()
# Purpose:
#   Run K independent replications of the full pipeline (46)+(47) → (50)+(53)+(54).
# Inputs:
#   r, eta, sigma : true parameters used to simulate each trajectory
#   T, dt         : time horizon and step size (N = round(T/dt))
#   K             : number of replications (default 50)
#   seeds         : optional integer vector of length K; if provided,
#                   replication k uses seeds[k] for reproducibility.
# Output:
#   A list with:
#     - results : K x 3 matrix with columns (sigma_hat, eta_hat, r_hat)
#     - means   : column means across K replications
#     - sds     : column standard deviations across K replications
#     - ci95    : 2-column matrix with mean ± 2*sd (a dispersion band)
#     - meta    : list with (r, eta, sigma, T, dt, K)
# Notes:
#   - Each replication draws Y0 from the Gamma law (47) unless y0 is forced
#     in run_replication (we keep the default: stationary init).
#   - The “±2·sd” bands summarize the spread of the estimator across replications;
#     they are NOT a confidence interval for the mean of the estimator.
# -----------------------------------------------------------------------------
run_montecarlo <- function(r, eta, sigma, T, dt, K = 50, seeds = NULL) {
  # storage for K triplets
  results <- matrix(NA_real_, nrow = K, ncol = 3L,
                    dimnames = list(NULL, c("sigma_hat", "eta_hat", "r_hat")))
  
  for (k in 1:K) {
    seed_k <- if (!is.null(seeds)) seeds[k] else k  # simple deterministic seeds if none provided
    rk <- run_replication(r = r, eta = eta, sigma = sigma,
                          T = T, dt = dt, seed = seed_k)
    results[k, ] <- unlist(rk)
  }
  
  # Monte Carlo summary statistics
  means <- colMeans(results)
  sds   <- apply(results, 2L, sd)
  ci95  <- cbind(lower = means - 2 * sds,
                 upper = means + 2 * sds)
  
  list(
    results = results,   # all K estimates (rows)
    means   = means,     # Monte Carlo averages
    sds     = sds,       # Monte Carlo standard deviations
    ci95    = ci95,      # mean ± 2*sd (dispersion band)
    meta    = list(r = r, eta = eta, sigma = sigma, T = T, dt = dt, K = K)
  )
}

# --------------------------
# Minimal usage example
# --------------------------
# This runs K=50 replications with modest T; increase T and/or K for tighter bands.
mc_out <- run_montecarlo(r = 2.0, eta = 0.5, sigma = 0.6,
                         T = 10, dt = 0.01, K = 50)

mc_out$means   # Monte Carlo means of (hat_sigma, hat_eta, hat_r)
mc_out$sds     # Monte Carlo standard deviations
mc_out$ci95    # Bands: mean ± 2*sd (spread of the estimator across replications)



# =========================================================
# Step 3 — Larger Monte Carlo scenarios (clear, English, equation-cited)
#
# Purpose:
#   Assess how estimator dispersion shrinks as the time horizon T and number
#   of replications K increase. Each scenario uses the pipeline:
#     (46)+(47)  simulation  →  (50) hat_sigma  →  (53) hat_eta  →  (54) hat_r
#   We report Monte Carlo mean, sd, and the band mean ± 2·sd for each estimator.
# =========================================================

# ---------- True parameters used to SIMULATE (46) with init (47) ----------
r_true     <- 2.0
eta_true   <- 0.5
sigma_true <- 0.6
true_vec   <- c(sigma_hat = sigma_true, eta_hat = eta_true, r_hat = r_true)

# ---------- Scenario 1 (larger than the tiny test): T=50, dt=0.01, K=200 ----------
T_big  <- 50
dt_big <- 0.01   # smaller dt (e.g., 0.005) reduces discretization bias but is slower
K_big  <- 200

mc_big <- run_montecarlo(r = r_true, eta = eta_true, sigma = sigma_true,
                         T = T_big, dt = dt_big, K = K_big)

# Summary vs. true values (means, sd, mean ± 2·sd, bias, relative bias)
summary_tab <- cbind(
  true     = true_vec[names(mc_big$means)],
  mean     = mc_big$means,
  sd       = mc_big$sds,
  lower    = mc_big$ci95[, "lower"],   # mean - 2*sd
  upper    = mc_big$ci95[, "upper"],   # mean + 2*sd
  bias     = mc_big$means - true_vec,
  rel_bias = (mc_big$means - true_vec) / true_vec
)
print(round(summary_tab, 6))


# ---------- Scenario 2: T=100, dt=0.01, K=200 ----------
T_2  <- 100
dt_2 <- 0.01
K_2  <- 200

mc_2 <- run_montecarlo(r = r_true, eta = eta_true, sigma = sigma_true,
                       T = T_2, dt = dt_2, K = K_2)

summary_2 <- cbind(
  true     = true_vec[names(mc_2$means)],
  mean     = mc_2$means,
  sd       = mc_2$sds,
  lower    = mc_2$ci95[, "lower"],
  upper    = mc_2$ci95[, "upper"],
  bias     = mc_2$means - true_vec,
  rel_bias = (mc_2$means - true_vec) / true_vec
)
cat("\n=== Summary — Scenario 2 (T=100, dt=0.01, K=200) ===\n")
print(round(summary_2, 6))


# ---------- Scenario 3 (verification, much larger): T=150, dt=0.005, K=500 ----------
T_3  <- 150
dt_3 <- 0.005
K_3  <- 500

mc_3 <- run_montecarlo(r = r_true, eta = eta_true, sigma = sigma_true,
                       T = T_3, dt = dt_3, K = K_3)

summary_3 <- cbind(
  true     = true_vec[names(mc_3$means)],
  mean     = mc_3$means,
  sd       = mc_3$sds,
  lower    = mc_3$ci95[, "lower"],
  upper    = mc_3$ci95[, "upper"],
  bias     = mc_3$means - true_vec,
  rel_bias = (mc_3$means - true_vec) / true_vec
)
cat("\n=== Summary — Scenario 3 (T=150, dt=0.005, K=500) ===\n")
print(round(summary_3, 6))



# =========================================================
# Step 4 — Comparative table for 3 Monte Carlo scenarios
# Purpose:
#   Summarize, side-by-side, the sampling behavior of the estimators
#   (hat_sigma from (50), hat_eta from (53), hat_r from (54)) under three
#   different (T, dt, K) settings using the pipeline (46)+(47) → (50)+(53)+(54).
# Notes:
#   - Uses only base R. No external packages required.
#   - Rounds only numeric columns to avoid errors with character columns.
# =========================================================

# Base packages presence check (friendly guard; both are base R)
if (!requireNamespace("stats", quietly = TRUE)) stop("'stats' should be available.")
if (!requireNamespace("graphics", quietly = TRUE)) stop("'graphics' should be available.")

# ---------- True parameters used to SIMULATE (46) with init (47) ----------
r_true     <- 2.0
eta_true   <- 0.5
sigma_true <- 0.6
true_vec   <- c(sigma_hat = sigma_true, eta_hat = eta_true, r_hat = r_true)

# -----------------------------------------------------------------------------
# run_and_summarize()
# Purpose:
#   Run a Monte Carlo scenario with given (T, dt, K), then assemble a tidy
#   data.frame with mean, sd, mean ± 2·sd, bias and relative bias for each
#   estimator (hat_sigma, hat_eta, hat_r).
# Inputs:
#   T_val, dt_val, K_val : scenario controls (time horizon, step, replications)
#   label                : human-readable scenario name for the table
# Output:
#   data.frame with columns:
#     scenario, parameter, true, mean, sd, lower, upper, bias, rel_bias, T, dt, K
# -----------------------------------------------------------------------------
run_and_summarize <- function(T_val, dt_val, K_val, label) {
  mc <- run_montecarlo(r = r_true, eta = eta_true, sigma = sigma_true,
                       T = T_val, dt = dt_val, K = K_val)
  
  df <- data.frame(
    scenario  = as.character(label),
    parameter = as.character(names(mc$means)),               # "sigma_hat","eta_hat","r_hat"
    true      = as.numeric(true_vec[names(mc$means)]),       # ground truth for each estimator
    mean      = as.numeric(mc$means),                        # MC mean
    sd        = as.numeric(mc$sds),                          # MC sd
    lower     = as.numeric(mc$ci95[, "lower"]),              # mean - 2*sd
    upper     = as.numeric(mc$ci95[, "upper"]),              # mean + 2*sd
    bias      = as.numeric(mc$means - true_vec[names(mc$means)]),
    rel_bias  = as.numeric((mc$means - true_vec[names(mc$means)]) / true_vec[names(mc$means)]),
    T         = as.numeric(T_val),
    dt        = as.numeric(dt_val),
    K         = as.integer(K_val),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  df
}

# ---------- Define the three scenarios ----------
tab1 <- run_and_summarize(T_val = 10,  dt_val = 0.01, K_val = 50,
                          label = "Scenario 1: T=10,  K=50,  dt=0.01")
tab2 <- run_and_summarize(T_val = 100, dt_val = 0.01, K_val = 200,
                          label = "Scenario 2: T=100, K=200, dt=0.01")
tab3 <- run_and_summarize(T_val = 150, dt_val = 0.005, K_val = 500,
                          label = "Scenario 3: T=150, K=500, dt=0.005")

tabla_final <- rbind(tab1, tab2, tab3)

# ---------- Pretty print: round numeric columns only ----------
nums <- sapply(tabla_final, is.numeric)
tabla_final_round <- tabla_final
tabla_final_round[, nums] <- round(tabla_final_round[, nums], 6)

print(tabla_final_round)


# -----------------------------------------------------------------------------
# plot_line_param()
# Purpose:
#   Draws the evolution of a given estimator (sigma_hat, eta_hat, r_hat)
#   across the Monte Carlo scenarios.
#
# Elements:
#   - X-axis: time horizon T (scenarios ordered by T).
#   - Black points/line: Monte Carlo mean of the estimator.
#   - Gray whiskers: mean ± 2*sd (sd = standard deviation across replications).
#   - Red dashed line: true parameter value.
#   - Numbers above points: scenario indices (1, 2, 3).
#   - Bottom legend: maps index → full scenario description (T, K, dt).
#   - Top-right legend: explains styles ("mean", "mean ± 2sd", "true value").
#
# Reading:
#   - Bias: closeness of black points to the red line.
#   - Variability: whisker length (shorter = more stable).
#   - Larger T (and/or smaller dt, larger K) should reduce whisker length and bias.
# -----------------------------------------------------------------------------

plot_line_param <- function(param_name, summary_table, true_vec) {
  df <- subset(summary_table, parameter == param_name)
  # order by T for left-to-right progression
  ord <- order(df$T); df <- df[ord, ]
  
  x <- df$T
  means  <- df$mean
  lowers <- df$lower
  uppers <- df$upper
  truep  <- true_vec[param_name]
  
  # Index labels (1,2,3) and mapping to scenarios
  idx_labels  <- paste0(1:nrow(df))
  scen_labels <- paste0(idx_labels, ": ", df$scenario)
  
  # Y-axis range with padding
  y_rng <- range(c(lowers, uppers, truep))
  y_pad <- 0.05 * diff(y_rng)
  y_rng[1] <- y_rng[1] - y_pad
  
  # Base plot
  plot(x, means, type = "o", pch = 19, ylim = y_rng,
       xlab = "Time horizon T", ylab = param_name,
       main = paste0("Evolution of ", param_name))
  
  # Whiskers = mean ± 2*sd
  arrows(x0 = x, y0 = lowers, x1 = x, y1 = uppers,
         angle = 90, code = 3, length = 0.06, col = "darkgray")
  
  # True value line
  abline(h = truep, lty = 2, lwd = 2, col = "red")
  
  # Numbers next to points
  text(x, means, labels = idx_labels, pos = 3, cex = 0.9)
  
  # Top-right legend (styles)
  legend("topright",
         legend = c("mean", "mean ± 2sd", "true value"),
         lty = c(1, 1, 2), pch = c(19, NA, NA),
         col = c("black", "darkgray", "red"),
         lwd = c(1, 1, 2), bty = "n", cex = 0.8)
  
  # Bottom legend (scenario mapping)
  x_mid <- mean(range(x))
  y_bot <- y_rng[1] + 0.01 * diff(y_rng)
  legend(x = x_mid, y = y_bot,
         legend = scen_labels,
         xjust = 0.5, yjust = 0, bty = "n", cex = 0.75)
}

# -----------------------------------------------------------------------------
# Draw the 3-panel figure (sigma_hat, eta_hat, r_hat side by side)
# -----------------------------------------------------------------------------
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3), mar = c(8, 4, 3, 1))
for (p in c("sigma_hat", "eta_hat", "r_hat")) {
  plot_line_param(p, tabla_final_round, true_vec)
}
par(op)


# -----------------------------------------------------------------------------
# plot_line_param_dual()
# Purpose:
#   Same as plot_line_param, but shows two types of bands:
#   (1) ±2*sd (gray)   = variability of the estimator across replications.
#   (2) ±2*sd/sqrt(K) (blue) = approximate 95% confidence interval for the
#                              Monte Carlo mean of the estimator.
# -----------------------------------------------------------------------------

plot_line_param_dual <- function(param_name, summary_table, true_vec) {
  df <- subset(summary_table, parameter == param_name)
  ord <- order(df$T); df <- df[ord, ]
  
  x <- df$T
  means  <- df$mean
  lowers <- df$lower
  uppers <- df$upper
  truep  <- true_vec[param_name]
  
  # K varies by scenario
  K_vals <- df$K
  
  # CI for the Monte Carlo mean (±2*sd/sqrt(K))
  lower_mean <- means - 2 * df$sd / sqrt(K_vals)
  upper_mean <- means + 2 * df$sd / sqrt(K_vals)
  
  # Scenario index labels
  idx_labels  <- paste0(1:nrow(df))
  scen_labels <- paste0(idx_labels, ": ", df$scenario)
  
  # Y-axis range with padding
  y_rng <- range(c(lowers, uppers, lower_mean, upper_mean, truep))
  y_pad <- 0.05 * diff(y_rng)
  y_rng[1] <- y_rng[1] - y_pad
  
  # Base plot
  plot(x, means, type = "o", pch = 19, ylim = y_rng,
       xlab = "Time horizon T", ylab = param_name,
       main = paste0("Evolution of ", param_name))
  
  # (1) ±2*sd whiskers (gray)
  arrows(x0 = x, y0 = lowers, x1 = x, y1 = uppers,
         angle = 90, code = 3, length = 0.06, col = "darkgray")
  
  # (2) ±2*sd/sqrt(K) whiskers (blue)
  arrows(x0 = x, y0 = lower_mean, x1 = x, y1 = upper_mean,
         angle = 90, code = 3, length = 0.06, col = "blue")
  
  # True value
  abline(h = truep, lty = 2, lwd = 2, col = "red")
  
  # Numbers above points
  text(x, means, labels = idx_labels, pos = 3, cex = 0.9)
  
  # Top-right legend
  legend("topright",
         legend = c("mean", "±2*sd (dispersion)", "±2*sd/sqrt(K) (MC mean CI)", "true value"),
         lty = c(1, 1, 1, 2), pch = c(19, NA, NA, NA),
         col = c("black", "darkgray", "blue", "red"),
         lwd = c(1, 1, 1, 2), bty = "n", cex = 0.75)
  
  # Bottom legend (scenario mapping)
  x_mid <- mean(range(x))
  y_bot <- y_rng[1] + 0.01 * diff(y_rng)
  legend(x = x_mid, y = y_bot,
         legend = scen_labels,
         xjust = 0.5, yjust = 0, bty = "n", cex = 0.75)
}

# Draw the 3-panel figure with dual bands
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3), mar = c(8, 4, 3, 1))
for (p in c("sigma_hat", "eta_hat", "r_hat")) {
  plot_line_param_dual(p, tabla_final_round, true_vec)
}
par(op)
