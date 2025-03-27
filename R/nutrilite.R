#' @title Estimate Nutrient Uptake Rates in Stream Ecosystems
#'
#' @description 
#' Nutrilite is a lightweight tool for estimating nutrient uptake rates in stream ecosystems 
#' using tracer injection experiments. It solves differential equations that model nutrient 
#' dynamics (advection, dispersion, and uptake) and fits these models to observed breakthrough 
#' curves (BTCs) of nutrient concentrations over time. By employing Bayesian optimization 
#' via \code{ParBayesianOptimization::bayesOpt}, Nutrilite efficiently infers globally optimal 
#' parameters for robust uptake rate estimation.
#' 
#' The function supports four model types:
#' \itemize{
#'   \item{\code{first_order}: Assumes uptake rate is proportional to nutrient concentration.}
#'   \item{\code{first_order_transit_storage}: Adds transient storage to the first-order model.}
#'   \item{\code{michaelis_menten}: Uses a Michaelis-Menten kinetics model for uptake saturation.}
#'   \item{\code{michaelis_menten_transit_storage}: Combines Michaelis-Menten kinetics with transient storage.}
#' }
#' 
#' This tool is particularly useful for hydrologists and ecologists studying nutrient cycling 
#' in streams, such as chloride (Cl) and nitrogen (N) dynamics.
#'
#' @param btc_data A data frame containing observed breakthrough curve data with the following columns:
#'   \itemize{
#'     \item{\code{cl}: Chloride concentration in mg/L.}
#'     \item{\code{n}: Nitrogen concentration in µg/L.}
#'     \item{\code{time}: Time since injection in seconds.}
#'   }
#' @param distance Numeric. Distance from the injection point to the sampling point, in meters.
#' @param cl_mass Numeric. Total mass of chloride injected, in milligrams.
#' @param n_mass Numeric. Total mass of nitrogen injected, in micrograms.
#' @param cross_area Numeric. Cross-sectional area of the stream at the sampling point, in square meters.
#' @param bounds Optional. A named list specifying lower and upper bounds for model parameters. 
#'   If \code{NULL} (default), bounds are automatically estimated from the data. Example structure:
#'   \preformatted{
#'   list(
#'     Cl = c(0, 10),      # Background Cl concentration (mg/L)
#'     N = c(0, 5),        # Background N concentration (µg/L)
#'     D = c(0.1, 10),     # Dispersion coefficient (m²/s)
#'     U = c(0.01, 2),     # Stream velocity (m/s)
#'     K = c(0.001, 0.1),  # First-order uptake rate constant (1/s)
#'     Alpha = c(0, 1),    # Transient storage exchange rate (1/s)
#'     AsA = c(0, 10),     # Ratio of storage zone to main channel area
#'     Vmax = c(0, 10000), # Maximum uptake rate for Michaelis-Menten (µg/L/s)
#'     Km = c(0, 1000000)  # Half-saturation constant for Michaelis-Menten (µg/L)
#'   )
#'   }
#' @param type Character string specifying the model type. One of:
#'   \code{"first_order"}, \code{"first_order_transit_storage"}, 
#'   \code{"michaelis_menten"}, or \code{"michaelis_menten_transit_storage"}. 
#'   Defaults to \code{"first_order"}.
#' @param initPoints Integer. Number of initial random evaluations for Bayesian optimization. 
#'   Default: \code{10}.
#' @param iters.n Integer. Number of Bayesian optimization iterations after initial sampling. 
#'   Default: \code{50}.
#' @param ... Additional arguments passed to \code{ParBayesianOptimization::bayesOpt}, 
#'   such as \code{acq} (acquisition function type: \code{"ei"}, \code{"ucb"}, \code{"poi"}), 
#'   \code{parallel} (logical for parallel evaluation), \code{kappa} (exploration-exploitation trade-off), 
#'   or \code{verbose} (output verbosity level).
#'
#' @return 
#' An object of class \code{bayesOpt} from the \code{ParBayesianOptimization} package, 
#' containing the results of the Bayesian optimization process. Key components include:
#' \itemize{
#'   \item{\code{FUN}: The scoring function used for optimization.}
#'   \item{\code{bounds}: The parameter bounds provided by the user.}
#'   \item{\code{iters}: Total number of iterations run.}
#'   \item{\code{scoreSummary}: A \code{data.table} with columns for parameter values, 
#'     scores (negative residual sum of squares), iteration details, and additional outputs from the scoring function.}
#'   \item{\code{stopStatus}: Reason for optimization termination (e.g., time limit reached, maximum iterations).}
#'   \item{\code{elapsedTime}: Total optimization time in seconds.}
#' }
#' To extract the optimal parameters, use:
#' \preformatted{
#' best_params <- result$scoreSummary[which.max(result$scoreSummary$Score), ]
#' }
#' For more details, see the \href{https://cran.r-project.org/web/packages/ParBayesianOptimization/index.html}{ParBayesianOptimization documentation}.
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(ParBayesianOptimization)
#' library(ReacTran)
#' 
#' # Simulate a realistic breakthrough curve dataset
#' set.seed(123)
#' time <- seq(0, 100, by = 5)  # Time in seconds
#' cl <- 0.1 + 10 * exp(-((time - 20)^2) / 50)  # Chloride concentration (mg/L)
#' n <- 2 + 50 * exp(-((time - 25)^2) / 60)     # Nitrogen concentration (µg/L)
#' btc_data <- data.frame(cl = cl, n = n, time = time)
#' 
#' # Set experimental parameters
#' distance <- 50      # Distance in meters
#' cl_mass <- 1000     # Chloride mass in mg
#' n_mass <- 500       # Nitrogen mass in µg
#' cross_area <- 0.5   # Cross-sectional area in m²
#' 
#' # Define parameter bounds (narrowed for faster optimization)
#' bounds <- list(
#'   Cl = c(0, 1),      # Background Cl concentration (mg/L)
#'   N = c(0, 5),       # Background N concentration (µg/L)
#'   D = c(0.1, 1),     # Dispersion coefficient (m²/s)
#'   U = c(0.1, 1),     # Stream velocity (m/s)
#'   K = c(0.001, 0.1)  # First-order uptake rate constant (1/s)
#' )
#' 
#' # Run with first-order model (limited iterations for quick execution)
#' result_fo <- nutrilite(btc_data, distance, cl_mass, n_mass, cross_area, 
#'                        bounds = bounds, type = "first_order", 
#'                        initPoints = 5, iters.n = 10)
#' 
#' # Extract best parameters from scoreSummary
#' best_params_fo <- result_fo$scoreSummary[which.max(result_fo$scoreSummary$Score), ]
#' print(best_params_fo)
#' }
#'
#' @export
#' @importFrom ReacTran setup.grid.1D fiadeiro tran.1D
#' @importFrom deSolve ode.1D
#' @importFrom ParBayesianOptimization bayesOpt

nutrilite <- function(
    btc_data, # concentration data with columns cl, n, time
    distance, # release distance 
    cl_mass, # injected Cl mass mg
    n_mass, # injected N mass ug
    cross_area, # cross_section_area m2
    bounds = NULL, # lower and upper bounds for parameters
    type = c('first_order', 'first_order_transit_storage', 'michaelis_menten', 'michaelis_menten_transit_storage'),
    initPoints = 10, # number of initial points
    iters.n = 50, # number of iterations
    ...) {

    # check if btc_data contains the required columns
    if(!all(c("cl", "n", "time") %in% colnames(btc_data))) {
        stop("btc_data must contain columns 'cl', 'n', and 'time'")
    }
    
    N_grid <-  ceiling((2 * distance) / 0.05) # number of grid cells
    grid <- ReacTran::setup.grid.1D(N = N_grid, dx.1 = 0.05)
    release_offset <- ceiling(N_grid / 40) # release offset ensures that the release point is not at the boundary
    release_zone <- ceiling(N_grid / 200) # release length
    injected_volume <- cross_area * release_zone * 0.05
    cl_injected <- cl_mass / injected_volume
    n_injected <- n_mass / injected_volume
    index <-  which.min(abs(grid$x.mid - distance)) + release_offset + 1 # index of the sampling point # the first col is time.

    # Advection-dispersion with first order uptake in main channel
    nut_up = function (times, y, parm){
        cs = y
        afdw = ReacTran::fiadeiro(v=parm[3], D=parm[2], grid=grid)
        trans = ReacTran::tran.1D(C=cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=afdw)
        uptake = parm[4] * cs
        input = parm[4] * parm[1]
        dcs = trans$dC - uptake + input
        list(c(dcs))
    }

    btc = function(parm){
        parm_cl = c(parm[[1]], parm[[3]], parm[[4]], 0)
        parm_n = c(parm[[2]], parm[[3]], parm[[4]], parm[[5]])
        yini_cl = c(rep(parm_cl[1], release_offset), rep(cl_injected, release_zone), rep(parm_cl[1], (N_grid - release_offset - release_zone)))
        yini_n = c(rep(parm_n[1], release_offset), rep(n_injected, release_zone), rep(parm_n[1], (N_grid - release_offset - release_zone)))
        cl_conc = deSolve::ode.1D(func=nut_up, y=yini_cl, parms=parm_cl, times=btc_data$time, nspec=1)[, index]
        n_conc = deSolve::ode.1D(func=nut_up, y=yini_n, parms=parm_n, times=btc_data$time, nspec=1)[, index]
        return(c(cl_conc, n_conc))
    }

    residual = function(Cl, N, D, U, K){
        parm = list(Cl, N, D, U, K)
        simu = btc(parm)
        resid = c(btc_data$cl, btc_data$n) - simu
        rss = sum(resid^2)
        return(list(Score = -rss))
    }

    # Advection-dispersion with transient storage and first order uptake in main channel
    nut_up_ts = function (times, y, parm){
        cs = y[1:N_grid]
        cts = y[(N_grid + 1):(2 * N_grid)]
        afdw = ReacTran::fiadeiro(v=parm[3], D=parm[2], grid=grid)
        trans = ReacTran::tran.1D(C=cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=afdw)
        ts_exchange = parm[5] * (cs - cts)
        uptake = parm[4] * cs
        input = parm[4] * parm[1]
        dcs = trans$dC - ts_exchange - uptake + input
        dcts = parm[6] * ts_exchange
        list(c(dcs, dcts))
    }

    # Fit Cl and N together, assuming same residual variance
    btc_ts = function(parm){
        parm_cl = c(parm[[1]], parm[[3]], parm[[4]], 0, parm[[5]], parm[[6]])
        parm_n = c(parm[[2]], parm[[3]], parm[[4]], parm[[7]], parm[[5]], parm[[6]])
        yini_cl = c(rep(parm_cl[1], release_offset), rep(cl_injected, release_zone), rep(parm_cl[1], (N_grid - release_offset - release_zone)), rep(parm_cl[1], N_grid))
        yini_n = c(rep(parm_n[1], release_offset), rep(n_injected, release_zone), rep(parm_n[1], (N_grid - release_offset - release_zone)), rep(parm_n[1], N_grid))
        cl_conc = deSolve::ode.1D(func=nut_up_ts, y=yini_cl, parms=parm_cl, times=btc_data$time, nspec=2)[, index]
        n_conc = deSolve::ode.1D(func=nut_up_ts, y=yini_n, parms=parm_n, times=btc_data$time, nspec=2)[, index]
        return(c(cl_conc, n_conc))
    }

    residual_ts = function(Cl, N, D, U, Alpha, AsA, K){
        parm = list(Cl, N, D, U, Alpha, AsA, K)
        simu = btc_ts(parm)
        resid = c(btc_data$cl, btc_data$n) - simu
        rss = sum(resid^2)
        return(list(Score = -rss))
    }

    # Advection-dispersion with M-M uptake in main channel
    nut_up_mm = function (times, y, parm){
        cs = y
        afdw = ReacTran::fiadeiro(v=parm[3], D=parm[2], grid=grid)
        trans = ReacTran::tran.1D(C=cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=afdw)
        uptake = parm[4] * cs / (cs + parm[5])
        input = parm[4] * parm[1] / (parm[1] + parm[5])
        dcs = trans$dC - uptake + input
        list(c(dcs))
    }

    btc_mm = function(parm){
        parm_cl = c(parm[[1]], parm[[3]], parm[[4]], 0, 0)
        parm_n = c(parm[[2]], parm[[3]], parm[[4]], parm[[5]], parm[[6]])
        yini_cl = c(rep(parm_cl[1], release_offset), rep(cl_injected, release_zone), rep(parm_cl[1], (N_grid - release_offset - release_zone))) 
        yini_n = c(rep(parm_n[1], release_offset), rep(n_injected, release_zone), rep(parm_n[1], (N_grid - release_offset - release_zone))) 
        cl_conc = deSolve::ode.1D(func=nut_up_mm, y=yini_cl, parms=parm_cl, times=btc_data$time, nspec=1)[, index]
        n_conc = deSolve::ode.1D(func=nut_up_mm, y=yini_n, parms=parm_n, times=btc_data$time, nspec=1)[, index]
        return(c(cl_conc, n_conc))
    }

    residual_mm = function(Cl, N, D, U, Vmax, Km){
        parm = list(Cl, N, D, U, Vmax, Km)
        simu = btc_mm(parm)
        resid = c(btc_data$cl, btc_data$n) - simu
        rss = sum(resid^2)
        return(list(Score = -rss))
    }

    # Advection-dispersion with transient storage and M-M uptake in main channel
    nut_up_ts_mm = function (times, y, parm){ 
        cs = y[1:N_grid]
        cts = y[(N_grid + 1):(2 * N_grid)]
        afdw = ReacTran::fiadeiro(v=parm[3], D=parm[2], grid=grid)
        trans = ReacTran::tran.1D(C=cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=afdw)
        ts_exchange = parm[5] * (cs - cts)
        uptake = parm[4] * cs / (cs + parm[7])
        input = parm[4] * parm[1] / (parm[1] + parm[7])
        dcs = trans$dC - ts_exchange - uptake + input
        dcts = parm[6] * ts_exchange
        list(c(dcs, dcts))
    }

    btc_ts_mm = function(parm){
        parm_cl = c(parm[[1]], parm[[3]], parm[[4]], 0, parm[[5]], parm[[6]], 0)
        parm_n = c(parm[[2]], parm[[3]], parm[[4]], parm[[7]], parm[[5]], parm[[6]], parm[[8]])
        yini_cl = c(rep(parm_cl[1], release_offset), rep(cl_injected, release_zone), rep(parm_cl[1], (N_grid - release_offset - release_zone)), rep(parm_cl[1], N_grid)) 
        yini_n = c(rep(parm_n[1], release_offset), rep(n_injected, release_zone), rep(parm_n[1], (N_grid - release_offset - release_zone)), rep(parm_n[1], N_grid)) 
        cl_conc = deSolve::ode.1D(func=nut_up_ts_mm, y=yini_cl, parms=parm_cl, times=btc_data$time, nspec=2)[, index]
        n_conc = deSolve::ode.1D(func=nut_up_ts_mm, y=yini_n, parms=parm_n, times=btc_data$time, nspec=2)[, index]
        return(c(cl_conc, n_conc))
    }

    residual_ts_mm = function(Cl, N, D, U, Alpha, AsA, Vmax, Km){
        parm = list(Cl, N, D, U, Alpha, AsA, Vmax, Km)
        simu = btc_ts_mm(parm)
        resid = c(btc_data$cl, btc_data$n) - simu
        rss = sum(resid^2)
        return(list(Score = -rss))
    }

    # if bounds are not provided, set them based on the estimated values
    if (is.null(bounds)) {
        # estimate initial parameters bounds
        # select the first 5 smallest values of cl and n as background
        n_bg <- mean(sort(btc_data$n[btc_data$n > 0])[1:5])
        cl_bg <- mean(sort(btc_data$cl[btc_data$cl > 0])[1:5])

        # Estimate the stream velocity U
        u_est <- distance / (btc_data$time[which.max(btc_data$cl)])
        
        # Estimate the dispersion coefficient D
        peak_cl <- max(btc_data$cl)
        peak_time_cl <- btc_data$time[which.max(btc_data$cl)]
        d_est <- (cl_mass^2) / (4 * pi * (peak_cl - cl_bg)^2 * cross_area^2 * peak_time_cl)

        # Estimate the first order rate constant K
        peak_n <- max(btc_data$n)
        peak_time_n <- btc_data$time[which.max(btc_data$n)]
        k_est <- -log((peak_n - n_bg) * (cross_area * sqrt(4 * pi * d_est * peak_time_n)) / n_mass) / peak_time_n

        lower_params <- list(
            Cl = cl_bg * 0.7,
            N = n_bg * 0.7,
            D = d_est * 0.5,
            U = u_est * 0.9,
            K = k_est * 0.7,
            Alpha = 0,
            AsA = 0,
            Vmax = 0,
            Km = 0
        )
        
        upper_params <- list(
            Cl = cl_bg * 1.3,
            N = n_bg * 1.3,
            D = d_est * 1.5,
            U = u_est * 1.1,
            K = k_est * 1.3,
            Alpha = 1,
            AsA = 10,
            Vmax = 10000,
            Km = 1000000
        )
        
        bounds <- list(
            Cl = c(lower_params$Cl, upper_params$Cl),
            N = c(lower_params$N, upper_params$N),
            D = c(lower_params$D, upper_params$D),
            U = c(lower_params$U, upper_params$U),
            K = c(lower_params$K, upper_params$K),
            Alpha = c(lower_params$Alpha, upper_params$Alpha),
            AsA = c(lower_params$AsA, upper_params$AsA),
            Vmax = c(lower_params$Vmax, upper_params$Vmax),
            Km = c(lower_params$Km, upper_params$Km)
        )
    }

    
    # choose one model to fit using package ParBayesianOptimization
    type <- tryCatch(
        match.arg(type),
        stop(
            "Invalid model type specified. Please choose one of the following: 'first_order', 'first_order_transit_storage', 'michaelis_menten', 'michaelis_menten_transit_storage'"
        )
    )
    opt_params <- switch(
        type,
        "first_order" = {
            ParBayesianOptimization::bayesOpt(
                Fun = residual,
                bounds = bounds[c('Cl', 'N', 'D', 'U', 'K')],
                initPoints = initPoints,
                iters.n = iters.n,
                ...
            )
        },
        "first_order_transit_storage" = {
            ParBayesianOptimization::bayesOpt(
                Fun = residual_ts,
                bounds = bounds[c('Cl', 'N', 'D', 'U', 'Alpha', 'AsA', 'K')],
                initPoints = initPoints,
                iters.n = iters.n,
                ...
            )
        },
        "michaelis_menten" = {
            ParBayesianOptimization::bayesOpt(
                Fun = residual_mm,
                bounds = bounds[c('Cl', 'N', 'D', 'U', 'Vmax', 'Km')],
                initPoints = initPoints,
                iters.n = iters.n,
                ...
            )
        },
        "michaelis_menten_transit_storage" = {
            ParBayesianOptimization::bayesOpt(
                Fun = residual_ts_mm,
                bounds = bounds[c('Cl', 'N', 'D', 'U', 'Alpha', 'AsA', 'Vmax', 'Km')],
                initPoints = initPoints,
                iters.n = iters.n,
                ...
            )
        }
    )

    return(opt_params)
}