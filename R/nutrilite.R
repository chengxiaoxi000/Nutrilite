#' @title Estimate Parameters in Stream Ecosystems
#'
#' @description 
#' Nutrilite fits stream solute dynamics models using Bayesian optimization via mlr3mbo.
#' It supports a two-step fitting process: first fitting hydraulic parameters using a 
#' conservative tracer (e.g., Cl), and then fitting biogeochemical parameters using a 
#' reactive tracer (e.g., N).
#'
#' @param time Numeric vector of observation times in seconds.
#' @param conc Numeric vector of observed concentrations (mg/L or µg/L).
#' @param distance Numeric. Distance from the injection point to the sampling point, in meters.
#' @param injected_mass Numeric. Total mass of the tracer injected.
#' @param cross_area Numeric. Cross-sectional area of the stream, in square meters.
#' @param tracer_type Character. Either \code{"conservative"} (for conservative tracers like Cl) 
#'        or \code{"reactive"} (for reactive tracers like N). Defaults to \code{"conservative"}.
#' @param model_type Character string specifying the model type. One of:
#'   \code{"FO"}, \code{"FO_TS"}, \code{"MM"}, or \code{"MM_TS"}. Defaults to \code{"FO"}.
#' @param bounds A named list specifying lower and upper bounds for parameters.
#' @param fixed_params A named list of parameters to keep constant. This is highly useful 
#'        when fitting reactive tracers where hydraulic parameters are already known.
#' @param n_evals Integer. Number of optimization evaluations for mlr3mbo. Default: \code{80}.
#'
#' @return An object of class \code{nutrilite_fit} containing the best parameters, MSE, RSS, and metadata.
#' 
#' @export
#' @importFrom ReacTran setup.grid.1D fiadeiro tran.1D
#' @importFrom deSolve ode.1D
#' @importFrom paradox ps p_dbl
#' @importFrom bbotk ObjectiveRFunDt trm OptimInstanceBatchSingleCrit opt
#' @importFrom mlr3mbo srlrn default_gp
#' @importFrom data.table data.table
nutrilite <- function(
    time, 
    conc, 
    distance, 
    injected_mass, 
    cross_area, 
    tracer_type = c("conservative", "reactive"), 
    model_type = c("FO", "FO_TS", "MM", "MM_TS"),
    bounds = NULL,
    fixed_params = list(),
    n_evals = 80
) {
  
  tracer_type <- match.arg(tracer_type)
  model_type <- match.arg(model_type)
  is_cons <- (tracer_type == "conservative")
  
  # 1. 网格与物理空间设置
  dx <- 0.05
  release_offset_physical <- 10.0
  release_zone_physical <- 0.4
  release_offset_cells <- ceiling(release_offset_physical / dx)
  release_zone_cells <- max(1, ceiling(release_zone_physical / dx))
  
  x_injection_start_absolute <- release_offset_cells * dx
  x_observation_absolute <- x_injection_start_absolute + distance
  downstream_buffer_physical <- max(20, distance * 0.75)
  total_domain_length_physical <- x_observation_absolute + downstream_buffer_physical
  
  N_grid <- max(100, ceiling(total_domain_length_physical / dx))
  grid <- ReacTran::setup.grid.1D(N = N_grid, dx.1 = dx)
  
  idx_in_grid_x_mid <- which.min(abs(grid$x.mid - x_observation_absolute))
  observation_column_index <- idx_in_grid_x_mid + 1
  
  injected_volume <- cross_area * release_zone_physical
  ini_conc <- injected_mass / injected_volume

  # 2. 定义内部 PDE 求解器
  river_solver_pde <- function(times, y, parms) {
    cs <- y[1:N_grid]
    cts <- y[(N_grid + 1):(2 * N_grid)]
    
    bg_conc   <- parms[["bg_conc"]]
    D_val     <- parms[["D"]]
    U_val     <- parms[["U"]]
    Alpha_val <- parms[["Alpha"]]
    AsA_coeff <- parms[["AsA_coeff"]]
    m_type    <- parms[["model_type"]]
    is_c      <- parms[["is_cons"]]
    
    afdw <- ReacTran::fiadeiro(v = U_val, D = D_val, grid = grid)
    trans <- ReacTran::tran.1D(C = cs, C.up = bg_conc, C.down = bg_conc,
                               D = D_val, v = U_val, dx = grid, AFDW = afdw)
    ts_exchange <- Alpha_val * (cs - cts)
    
    uptake <- 0
    input_term <- 0
    if (!is_c) {
      if (m_type == "FO" || m_type == "FO_TS") {
        K_val <- parms[["K_N"]]
        uptake <- K_val * cs
        input_term <- K_val * bg_conc
      } else if (m_type == "MM" || m_type == "MM_TS") {
        Vmax_val <- parms[["V_max"]]
        Km_val   <- parms[["K_m"]]
        uptake <- Vmax_val * cs / (Km_val + cs)
        input_term <- Vmax_val * bg_conc / (Km_val + bg_conc)
      }
    }
    
    dcs <- trans$dC - ts_exchange - uptake + input_term
    dcts <- AsA_coeff * ts_exchange
    list(c(dcs, dcts))
  }

  # 3. 动态构建参数搜索空间
  param_list_dynamic <- list()
  if (is.null(bounds)) stop("Please provide a 'bounds' list.")
  
  if (is_cons) {
    param_list_dynamic$D <- paradox::p_dbl(lower = bounds$D[1], upper = bounds$D[2])
    param_list_dynamic$U <- paradox::p_dbl(lower = bounds$U[1], upper = bounds$U[2])
    param_list_dynamic$bg_conc <- paradox::p_dbl(lower = bounds$bg_conc[1], upper = bounds$bg_conc[2])
    if (grepl("TS", model_type)) {
      param_list_dynamic$Alpha <- paradox::p_dbl(lower = bounds$Alpha[1], upper = bounds$Alpha[2])
      param_list_dynamic$AsA_coeff <- paradox::p_dbl(lower = bounds$AsA_coeff[1], upper = bounds$AsA_coeff[2])
    }
  } else {
    param_list_dynamic$bg_conc <- paradox::p_dbl(lower = bounds$bg_conc[1], upper = bounds$bg_conc[2])
    if (grepl("FO", model_type)) {
      param_list_dynamic$K_N <- paradox::p_dbl(lower = bounds$K_N[1], upper = bounds$K_N[2])
    } else if (grepl("MM", model_type)) {
      param_list_dynamic$V_max <- paradox::p_dbl(lower = bounds$V_max[1], upper = bounds$V_max[2])
      param_list_dynamic$K_m <- paradox::p_dbl(lower = bounds$K_m[1], upper = bounds$K_m[2])
    }
  }

  # 4. 目标函数
  objective_fun <- function(xdt) {
    scores <- vapply(1:nrow(xdt), function(i) {
      params_optimizing_list <- as.list(xdt[i, ])
      current_full_params <- c(params_optimizing_list, fixed_params)
      
      possible_keys <- c("D", "U", "Alpha", "AsA_coeff", "bg_conc", "K_N", "V_max", "K_m")
      for(key in possible_keys) {
        if(is.null(current_full_params[[key]])) current_full_params[[key]] <- 0
        if(key == "AsA_coeff" && current_full_params[[key]] == 0) current_full_params[[key]] <- 1.0
      }
      
      parm_for_solver <- c(current_full_params, list(model_type = model_type, is_cons = is_cons))
      bg_val <- current_full_params$bg_conc
      
      yini <- c(rep(bg_val, release_offset_cells),
                rep(ini_conc, release_zone_cells),
                rep(bg_val, N_grid - release_offset_cells - release_zone_cells),
                rep(bg_val, N_grid))
                
      simulated_btc <- try(deSolve::ode.1D(func=river_solver_pde, y=yini, parms=parm_for_solver, times=time, nspec=2), silent = TRUE)
      current_mse <- 1e12
      
      if (!inherits(simulated_btc, "try-error")) {
        sim_vals <- simulated_btc[, observation_column_index + 1]
        if (length(sim_vals) == length(conc) && !any(is.na(sim_vals))) {
          rss_val <- sum((conc - sim_vals)^2)
          current_mse <- rss_val / length(conc)
        }
      }
      if (is.na(current_mse) || is.infinite(current_mse)) current_mse <- 1e11
      return(current_mse)
    }, FUN.VALUE = numeric(1))
    
    return(data.table::data.table(Score = scores))
  }

  # 5. 运行 mlr3mbo 优化
  best_params_optim <- list()
  best_mse <- NA
  
  if (length(param_list_dynamic) > 0) {
    param_set_dynamic <- do.call(paradox::ps, param_list_dynamic)
    objective_mlr3 <- bbotk::ObjectiveRFunDt$new(
      fun = objective_fun,
      domain = param_set_dynamic,
      codomain = paradox::ps(Score = paradox::p_dbl(tags = "minimize")),
      properties = "noisy"
    )
    terminator <- bbotk::trm("evals", n_evals = n_evals)
    instance <- bbotk::OptimInstanceBatchSingleCrit$new(objective = objective_mlr3, terminator = terminator)
    optimizer <- bbotk::opt("mbo", surrogate = mlr3mbo::srlrn(mlr3mbo::default_gp(noisy = TRUE)))
    
    optimizer$optimize(instance)
    
    best_params_optim <- as.list(instance$result_x_search_space)
    best_mse <- instance$result_y
  }

  final_all_params <- c(fixed_params, best_params_optim)
  rss_final <- if (!is.na(best_mse)) best_mse * length(conc) else NA

  # 6. 返回结果 (S3 Object)
  result <- list(
    best_params = final_all_params,
    best_mse = best_mse,
    rss = rss_final,
    n_data = length(conc),
    model_type = model_type,
    tracer_type = tracer_type,
    n_params_fitted = length(param_list_dynamic),
    data = data.frame(time = time, conc = conc),
    setup = list(distance = distance, injected_mass = injected_mass, cross_area = cross_area)
  )
  
  class(result) <- "nutrilite_fit"
  return(result)
}