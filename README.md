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
```R
# Run a minimal replication
res <- run_replication(r = 2.0, eta = 0.5, sigma = 0.6,
                       T = 10, dt = 0.01, seed = 1)
print(res)
