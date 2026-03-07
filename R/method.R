#' @export
predict.nutrilite_fit <- function(object, times = NULL, ...) {
  if (is.null(times)) {
    times <- object$data$time
  }
  
  sim_result <- simulate_btc(
    time_seq = times,
    distance = object$setup$distance,
    injected_mass = object$setup$injected_mass,
    cross_area = object$setup$cross_area,
    params = object$best_params,
    tracer_type = object$tracer_type,
    model_type = object$model_type
  )
  
  return(sim_result$conc)
}


#' @export
plot.nutrilite_fit <- function(x, ...) {
  obs_time <- x$data$time
  obs_conc <- x$data$conc
  
  # 生成200个平滑的时间点用于绘制预测曲线
  sim_time <- seq(min(obs_time), max(obs_time), length.out = 200)
  sim_conc <- predict(x, times = sim_time)
  
  tracer_label <- ifelse(x$tracer_type == "conservative", "(Conservative)", "(Reactive)")
  title <- paste("Nutrilite Fit:", x$model_type, tracer_label)
  
  # 基础 R 绘图
  plot(obs_time, obs_conc, type = "p", pch = 16, col = "black",
       xlab = "Time (seconds)", ylab = "Concentration",
       main = title)
  lines(sim_time, sim_conc, col = "red", lwd = 2)
  
  legend("topright", legend = c("Observed Data", "Fitted Model"),
         col = c("black", "red"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2))
}

#' @export
AIC.nutrilite_fit <- function(object, ..., k = 2) {
  n <- object$n_data
  rss <- object$rss
  n_params <- object$n_params_fitted
  
  if (is.na(rss) || is.na(n) || n == 0) {
    warning("Cannot calculate AIC: Missing RSS or data points.")
    return(NA_real_)
  }
  
  # 标准的最小二乘法 AIC 近似计算公式
  aic_value <- k * n_params + n * log(rss / n)
  
  return(aic_value)
}

#' @export
summary.nutrilite_fit <- function(object, ...) {
  cat("========================================\n")
  cat("          Nutrilite Model Fit           \n")
  cat("========================================\n")
  cat("Tracer Type   :", object$tracer_type, "\n")
  cat("Model Type    :", object$model_type, "\n")
  cat("Data Points   :", object$n_data, "\n")
  cat("Mean Sq Error :", signif(object$best_mse, 5), "\n")
  # 👇 这里直接调用原生的 AIC()
  cat("AIC           :", signif(AIC(object), 5), "\n") 
  cat("----------------------------------------\n")
  cat("Best Optimized Parameters:\n")
  
  param_df <- data.frame(Value = unlist(object$best_params))
  print(param_df)
  cat("========================================\n")
}