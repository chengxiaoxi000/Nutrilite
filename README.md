# Nutrilite: Estimate Nutrient Uptake in Stream Ecosystems

## Overview
**Nutrilite** is a lightweight R package designed for estimating nutrient uptake rates in stream ecosystems. It works by using numerical solutions of differential equations that describe nutrient dynamics, and then fitting these models to observed breakthrough curves (BTCs) of nutrient concentrations. 

By leveraging Bayesian optimization (via `mlr3mbo`), Nutrilite efficiently infers the globally optimal model parameters, enabling robust estimation of nutrient uptake rates.

## Key Features
* **Multiple Kinetic Models**: Supports First-Order (`FO`), First-Order with Transient Storage (`FO_TS`), Michaelis-Menten (`MM`), and Michaelis-Menten with Transient Storage (`MM_TS`) models.
* **Flexible Parameter Fitting**: Users can easily define which parameters to optimize (by providing `bounds`) and which to keep constant (by providing `fixed_params`).
* **BTC Simulator**: Includes a `simulate_btc()` function to generate simulated breakthrough curves. This is highly useful for experimental design (e.g., predicting peak arrival time and concentration) or sensitivity analysis before going to the field.
* **S3 Methods Integration**: Provides native S3 methods including `plot()`, `summary()`, `predict()`, and `AIC()` for easy model evaluation and visualization.

## Installation
Currently, the package is in development (version 0.0.0.9000). You can install the development version from GitHub using `devtools` or `remotes`:

```R
# install.packages("devtools")
# Replace "chengxiaoxi000" with your actual GitHub username if different
devtools::install_github("chengxiaoxi000/nutrilite")
```

## Quick Start

Nutrilite recommends a **two-step fitting workflow** to estimate the uptake rates of reactive nutrients (like Nitrogen or Phosphorus). The following example uses the built-in simulated dataset to demonstrate this process.

### Step 1: Fitting the Conservative Tracer
First, we fit the conservative tracer (e.g., Chloride) to estimate the stream's physical transport properties (such as Dispersion $D$, Velocity $U$, and Transient Storage exchange coefficient $\alpha$).

```R
library(Nutrilite)

# Load the built-in simulated conservative tracer (Cl) data
cl_file <- system.file("extdata", "simu_fo_ts_cl.csv", package = "Nutrilite")
df_cl <- read.csv(cl_file)

# Fit the conservative tracer
fit_cl <- nutrilite(
  time = df_cl$time, 
  conc = df_cl$cl_obs, 
  distance = 40,               # Distance from injection to sampling point in meters
  injected_mass = 1500 * 0.4,  # Injected mass
  cross_area = 1.0,            # Cross-sectional area of the stream
  tracer_type = "conservative", 
  model_type = "FO_TS", 
  bounds = list(
    D = c(0.1, 1),        
    U = c(1.35, 1.65),        
    Alpha = c(0.01, 0.1),   
    AsA_coeff = c(0.5, 3.5),
    bg_conc = c(5.5, 6.5)    
  ),
  n_evals = 80
)

# View the estimated physical and hydrological parameters
summary(fit_cl)
plot(fit_cl, main="Conservative Tracer Fit (Cl)")
```

### Step 2: Fitting the Reactive Tracer
Next, extract the optimized hydrological parameters from Step 1 and pass them as `fixed_params`. This isolates the reactive nutrient dynamics, allowing us to estimate the biogeochemical uptake rate (e.g., $K_N$).

```R
# Load the built-in simulated reactive tracer (N) data
n_file <- system.file("extdata", "simu_fo_ts_n.csv", package = "Nutrilite")
df_n <- read.csv(n_file)

# Extract physical parameters from the conservative fit
hydro_params <- fit_cl$best_params

# Fit the reactive tracer
fit_n <- nutrilite(
  time = df_n$time, 
  conc = df_n$n_obs, 
  distance = 40, 
  injected_mass = 2500 * 0.4, 
  cross_area = 1.0,
  tracer_type = "reactive", 
  model_type = "FO_TS", 
  # Fix the hydrological parameters estimated in Step 1
  fixed_params = list(
    D = hydro_params$D, 
    U = hydro_params$U, 
    Alpha = hydro_params$Alpha, 
    AsA_coeff = hydro_params$AsA_coeff
  ),
  # Only optimize the uptake rate (K_N) and background concentration
  bounds = list(
    K_N = c(0.01, 0.1),     
    bg_conc = c(2.5, 3.5)   
  ),
  n_evals = 80
)

# View the estimated uptake rate and model evaluation
summary(fit_n)
plot(fit_n, main="Reactive Tracer Fit (N)")
```

## Authors
* **Traci Chan** (chengx2023@lzu.edu.cn)
* **Chao Song** (chaosong@lzu.edu.cn)

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
