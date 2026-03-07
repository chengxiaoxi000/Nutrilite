#' @title Simulate a Breakthrough Curve (BTC)
#' @description 
#' Generates a simulated breakthrough curve based on user-defined parameters.
#' 
#' @param time_seq Numeric vector of times to simulate.
#' @param distance Distance in meters.
#' @param injected_mass Injected mass (mg or ug).
#' @param cross_area Cross-sectional area in square meters.
#' @param params A named list of model parameters (e.g., U, D, Alpha, AsA_coeff, bg_conc, K_N, V_max, K_m).
#' @param tracer_type "conservative" or "reactive".
#' @param model_type "FO", "FO_TS", "MM", or "MM_TS".
#' @return A data.frame with columns `time` and `conc`.
#' @export
simulate_btc <- function(time_seq, distance, injected_mass, cross_area, 
                         params, tracer_type = c("conservative", "reactive"), 
                         model_type = c("FO", "FO_TS", "MM", "MM_TS")) {
  
  tracer_type <- match.arg(tracer_type)
  model_type <- match.arg(model_type)
  is_cons <- (tracer_type == "conservative")
  
  # --- 1. 构建物理网格 (与 nutrilite 内部保持一致) ---
  dx <- 0.05
  release_offset_cells <- ceiling(10.0 / dx)
  release_zone_cells <- max(1, ceiling(0.4 / dx))
  
  x_injection_start <- release_offset_cells * dx
  x_observation <- x_injection_start + distance
  total_length <- x_observation + max(20, distance * 0.75)
  
  N_grid <- max(100, ceiling(total_length / dx))
  grid <- ReacTran::setup.grid.1D(N = N_grid, dx.1 = dx)
  obs_index <- which.min(abs(grid$x.mid - x_observation)) + 1
  
  injected_volume <- cross_area * 0.4
  ini_conc <- injected_mass / injected_volume
  
  # 提取参数 (缺失则补0或默认值)
  p <- list(
    bg_conc = ifelse(is.null(params$bg_conc), 0, params$bg_conc),
    U = ifelse(is.null(params$U), 0.1, params$U),
    D = ifelse(is.null(params$D), 0.1, params$D),
    Alpha = ifelse(is.null(params$Alpha), 0, params$Alpha),
    AsA_coeff = ifelse(is.null(params$AsA_coeff), 1.0, params$AsA_coeff),
    K_N = ifelse(is.null(params$K_N), 0, params$K_N),
    V_max = ifelse(is.null(params$V_max), 0, params$V_max),
    K_m = ifelse(is.null(params$K_m), 1, params$K_m),
    model_type = model_type,
    is_cons = is_cons
  )
  
  # --- 2. PDE 求解器 ---
  river_solver <- function(times, y, parms) {
    cs <- y[1:N_grid]
    cts <- y[(N_grid + 1):(2 * N_grid)]
    
    afdw <- ReacTran::fiadeiro(v = parms$U, D = parms$D, grid = grid)
    trans <- ReacTran::tran.1D(C = cs, C.up = parms$bg_conc, C.down = parms$bg_conc,
                               D = parms$D, v = parms$U, dx = grid, AFDW = afdw)
    ts_exchange <- parms$Alpha * (cs - cts)
    
    uptake <- 0; input_term <- 0
    if (!parms$is_cons) {
      if (grepl("FO", parms$model_type)) {
        uptake <- parms$K_N * cs
        input_term <- parms$K_N * parms$bg_conc
      } else if (grepl("MM", parms$model_type)) {
        uptake <- parms$V_max * cs / (parms$K_m + cs)
        input_term <- parms$V_max * parms$bg_conc / (parms$K_m + parms$bg_conc)
      }
    }
    
    dcs <- trans$dC - ts_exchange - uptake + input_term
    dcts <- parms$AsA_coeff * ts_exchange
    list(c(dcs, dcts))
  }
  
  # --- 3. 初始条件与求解 ---
  yini <- c(rep(p$bg_conc, release_offset_cells),
            rep(ini_conc, release_zone_cells),
            rep(p$bg_conc, N_grid - release_offset_cells - release_zone_cells),
            rep(p$bg_conc, N_grid))
            
  out <- deSolve::ode.1D(func = river_solver, y = yini, parms = p, times = time_seq, nspec = 2)
  sim_conc <- out[, obs_index + 1]
  
  return(data.frame(time = time_seq, conc = sim_conc))
}