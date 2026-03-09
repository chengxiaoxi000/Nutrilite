#' @title Simulate a Breakthrough Curve (BTC)
#' @description 
#' Generates a simulated breakthrough curve based on user-defined parameters.
#' This is highly useful for experimental design (e.g., predicting peak arrival time 
#' and concentration before going to the field) or sensitivity analysis.
#' 
#' @param time_seq Numeric vector of times (in seconds) to simulate.
#' @param distance Distance from injection to observation point in meters.
#' @param injected_mass Injected mass of the tracer (mg or ug).
#' @param cross_area Cross-sectional area of the stream in square meters.
#' @param params A named list of model parameters (e.g., U, D, Alpha, AsA_coeff, bg_conc, K_N, V_max, K_m).
#' @param tracer_type Character. Either "conservative" or "reactive". Defaults to "conservative".
#' @param model_type Character. "FO", "FO_TS", "MM", or "MM_TS". Defaults to "FO".
#' @return A data.frame with columns `time` and `conc`.
#' @export
#' @importFrom ReacTran setup.grid.1D fiadeiro tran.1D
#' @importFrom deSolve ode.1D
simulate_btc <- function(
    time_seq, 
    distance, 
    injected_mass, 
    cross_area, 
    params, 
    tracer_type = c("conservative", "reactive"), 
    model_type = c("FO", "FO_TS", "MM", "MM_TS")
) {
  
  tracer_type <- match.arg(tracer_type)
  model_type <- match.arg(model_type)
  is_cons <- (tracer_type == "conservative")
  
  # --- 1. 构建物理网格 (已与 nutrilite 主函数同步) ---
  dx <- 0.05
  
  # 动态设定上游缓冲带，保底10米，随距离略微增加
  release_offset_physical <- max(10.0, distance * 0.15)
  release_zone_physical <- 0.4
  
  release_offset_cells <- ceiling(release_offset_physical / dx)
  release_zone_cells <- max(1, ceiling(release_zone_physical / dx))
  
  x_injection_start_absolute <- release_offset_cells * dx
  x_observation_absolute <- x_injection_start_absolute + distance
  downstream_buffer_physical <- max(20, distance * 0.75)
  total_domain_length_physical <- x_observation_absolute + downstream_buffer_physical
  
  N_grid <- max(100, ceiling(total_domain_length_physical / dx))
  grid <- ReacTran::setup.grid.1D(N = N_grid, dx.1 = dx)
  
  obs_index <- which.min(abs(grid$x.mid - x_observation_absolute)) + 1
  
  injected_volume <- cross_area * release_zone_physical
  ini_conc <- injected_mass / injected_volume
  
  # --- 2. 提取并补齐参数 (赋予安全默认值) ---
  p <- list(
    bg_conc   = ifelse(is.null(params$bg_conc), 0, params$bg_conc),
    U         = ifelse(is.null(params$U), 0.1, params$U),
    D         = ifelse(is.null(params$D), 0.1, params$D),
    Alpha     = ifelse(is.null(params$Alpha), 0, params$Alpha),
    AsA_coeff = ifelse(is.null(params$AsA_coeff), 1.0, params$AsA_coeff),
    K_N       = ifelse(is.null(params$K_N), 0, params$K_N),
    V_max     = ifelse(is.null(params$V_max), 0, params$V_max),
    K_m       = ifelse(is.null(params$K_m), 1.0, params$K_m),
    model_type = model_type,
    is_cons    = is_cons
  )
  
  # --- 3. PDE 求解器 ---
  river_solver <- function(times, y, parms) {
    cs <- y[1:N_grid]
    cts <- y[(N_grid + 1):(2 * N_grid)]
    
    afdw <- ReacTran::fiadeiro(v = parms$U, D = parms$D, grid = grid)
    trans <- ReacTran::tran.1D(C = cs, C.up = parms$bg_conc, C.down = parms$bg_conc,
                               D = parms$D, v = parms$U, dx = grid, AFDW = afdw)
    ts_exchange <- parms$Alpha * (cs - cts)
    
    uptake <- 0
    input_term <- 0
    if (!parms$is_cons) {
      if (parms$model_type == "FO" || parms$model_type == "FO_TS") {
        uptake <- parms$K_N * cs
        input_term <- parms$K_N * parms$bg_conc
      } else if (parms$model_type == "MM" || parms$model_type == "MM_TS") {
        uptake <- parms$V_max * cs / (parms$K_m + cs)
        input_term <- parms$V_max * parms$bg_conc / (parms$K_m + parms$bg_conc)
      }
    }
    
    dcs <- trans$dC - ts_exchange - uptake + input_term
    dcts <- parms$AsA_coeff * ts_exchange
    list(c(dcs, dcts))
  }
  
  # --- 4. 初始条件与求解 ---
  yini <- c(rep(p$bg_conc, release_offset_cells),
            rep(ini_conc, release_zone_cells),
            rep(p$bg_conc, N_grid - release_offset_cells - release_zone_cells),
            rep(p$bg_conc, N_grid))
            
  out <- deSolve::ode.1D(func = river_solver, y = yini, parms = p, times = time_seq, nspec = 2)
  sim_conc <- out[, obs_index + 1]
  
  return(data.frame(time = time_seq, conc = sim_conc))
}