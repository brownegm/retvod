#' Retrieve VOD for all values
#'
#' @param tbH Horizontal brightness temperatures
#' @param tbV Vertical brightness temperatures
#' @param smc Soil moisture
#' @param vod VOD estimate
#' @param Tair Air temperature
#' @param Tsoil Soil temperature
#' @param omega Canopy opacity, scattering albedo
#' @param roughness Soil roughness estimate
#' @param inc_angle Incidence angle
#' @param cf Clay fraction
#' @param silent Silence progress bars in the console. Default is FALSE.
#'
#' @return Retrieved VOD and auxillary information
#' @export
#'
#' @import cli

retVOD <- function(tbH, tbV,
                   smc,
                   vod,
                   Tair,
                   Tsoil,
                   omega,
                   roughness,
                   cf,
                   inc_angle, silent = F) {
  if (length(tbH) != length(tbV) || length(tbH) != length(smc)) {
    stop("tbH, tbV, and smc lengths differ.")
  }

  ## calculate gamma for each VOD test value
  gamma <- exp(-vod / cos(inc_angle* (pi/180)))

  ## calculate epsilon (dielectric) and reflectivitys for each value of soil moisture
  eps_list <- sapply(smc, \(s) mironov(1.4e9, s, cf)$dielectric)
  reflecs <- sapply(eps_list, \(e) fresnelr(eps = e, theta = inc_angle, h=roughness), simplify = F)

  ## prepare omega vector (if only one value, repeat for each tbH)
  o <- if (length(omega) == 1) rep_len(omega, length(tbH)) else omega

  ## initialize lists for the results
  results <- vector("list", length(tbH))
  # names(results) <- c(
  #   "vodEst", "cfEst", "smEst",
  #   "tbHpred", "tbVpred",
  #   "tbHcost", "tbVcost", "res_mat"
  # )

  n <- length(tbH)

  results <- list(
    vodEst = numeric(n), cfEst = numeric(n), smEst = numeric(n),
    tbHpred = numeric(n), tbVpred = numeric(n),
    tbHcost = numeric(n), tbVcost = numeric(n), res_mat= numeric(n)
  )

  # for each brightness temperature retrieve VOD and soil moisture
  if (!silent) {
    cli::cli_progress_bar("Retrieving VOD values...", total = length(tbH))
  }

  for (i in seq_along(tbH)) {
    est <- solveSmVod(
      reflec = reflecs[i], tbH = tbH[i], tbV = tbV[i],
      gamma = gamma,# vod = vod,
      Tair = Tair[i], Tsoil = Tsoil[i],
      omega = o[i], mat = F
    )

    smidx<-which(sapply(reflecs[i], function(x) x$fH == est$reflec_best$fH && x$fV == est$reflec_best$fV))


    # return elements from retrieval
    results$vodEst[i] <- vod[which(gamma==est$gamma_best)]
    results$cfEst[i] <- est$cf_tb
    results$smEst[i] <- smc[i]#[smidx]
    results$tbVpred[i] <- est$pred_tbV
    results$tbHpred[i] <- est$pred_tbH
    results$tbHcost[i] <- est$cf_tbH
    results$tbVcost[i] <- est$cf_tbV
    results$res_mat[i] <- est$cf_mat

    if (!silent) {
      cli_progress_update() # Update progress bar
    }
  }

  results$rmse_k <- mean(sqrt(results$cfEst), na.rm = T)
  results$reflectivity <- reflecs
  results$gamma <- gamma

  if (!silent) {
    cli_progress_done() # Complete progress bar
    cli_alert_success("VOD retrieval complete!")
  }

  if (!identical(results$smEst, smc)) {
    warning("'Estimated' soil moisture values do not match input.")
  }

  return(structure(results,
    creation_time = Sys.time(),
    inputs = list(h = tbH, v = tbV, sm = smc, Tair = Tair, Tsoil=Tsoil, omega = omega, rough = roughness, angle = inc_angle),
    class = c( "retVOD","data.frame")
  ))
}


#' @export
#' @importFrom graphics hist par abline
#' @importFrom grDevices dev.new
plot.retVOD <- function(x, ...) {
  dev.new(width = 12, height = 8)

  inputs <- attr(x, "inputs")
  par(mfrow = c(3, 4), mar = c(5, 5, 1, 1))

  plot(x$vodEst,
    main = "Retrieved Vegetation Optical Depth",
    ylab = "VOD",
    xlab = "Index"
  )

  plot(x$cfEst,
    main = "Total Tb Residuals",
    ylab = "Brightness Temps (K^2)",
    xlab = "Index"
  )

  hist(x$cfEst, main = "Histogram Tb Residuals", xlab = "Tb Residuals")

  plot(x$tbHcost,
    main = "TbH Residuals (K^2)",
    ylab = "Brightness Temps (K^2)",
    xlab = "Index"
  )

  plot(x$tbVcost,
    main = "TbV Residuals (K^2)",
    ylab = "Brightness Temps (K^2)",
    xlab = "Index"
  )

  plot(x$tbHpred ~ x$tbHcost,
    main = "Predicted TbH ~ TbH Residuals",
    ylab = "Brightness Temps (K)",
    xlab = "Tb H Residuals (K^2)"
  )
  abline(a = 0, b = -1, col = "red")

  plot(x$tbVpred ~ x$tbVcost,
    main = "Predicted TbV ~ TbV Residuals",
    ylab = "Brightness Temps (K)",
    xlab = "Tb V Residuals (K^2)"
  )
  abline(a = 0, b = -1, col = "red")

  plot(inputs$h ~ x$tbHcost,
    main = "Observed TbH ~ TbH Residuals",
    ylab = "Brightness Temps (K)",
    xlab = "Tb H Residuals (K^2)"
  )
  abline(a = 0, b = 1, col = "red")
  plot(inputs$v ~ x$tbVcost,
    main = "Observed TbV ~ TbV Residuals",
    ylab = "Brightness Temps (K)",
    xlab = "Tb V Residuals (K^2)"
  )
  abline(a = 0, b = 1, col = "red")

  plot(x$cfEst ~ inputs$sm,
       main = "Total Res. ~ Soil Moisture",
       ylab = "Total Tb Residuals",
       xlab = "Soil Moisture")

  plot(x$cfEst ~ inputs$Tair,
       main = "Total Res. ~ Air Temp",
       ylab = "Total Tb Residuals",
       xlab = "Air Temperature")

  plot(x$cfEst ~ inputs$Tsoil,
       main = "Total Res. ~ Soil Temp",
       ylab = "Total Tb Residuals",
       xlab = "Soil Temperature")

}
