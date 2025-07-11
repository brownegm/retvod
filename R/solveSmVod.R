#' Solve for soil moisture and VOD
#'
#' @param reflec H and V pol reflectivities
#' @param gamma Gamma estimates for each test VOD value
#' @param tbH Horizontal brightness temperature
#' @param tbV Vertical brightness temperature
#' @param Tair Air temperature
#' @param Tsoil Soil temperature
#' @param omega Scattering albedo
#' @param tno TRUE if omega is to be varied, FALSE if it is constant
#' @param mat TRUE returns the cost function matrix
#'
#' @description This function uses the a range of input soil moisture and vod values to solve for the best VOD value given the air, soil, and observed brightness temperatures.
#' @details The function returns the predicted brightness temperatures (i.e., using tau-omega), residuals for each polarization, and the predicted soil moisture and VOD. Mironov is used to determine the dielectric constant for a given soil moisture and has a set frequency of 1.4e9Hz (for L-band). Also the clay fraction at MOFLUX is 23.2%.
#' @return List of predicted brightness temperatures, soil moisture and VOD:
#'
#' @export
#'
solveSmVod <- function(reflec,
                       gamma,
                       tbH, tbV,
                       Tair, Tsoil,
                       omega, mat=F, tno = F) {

  ## initialize output matrices
  grid <- expand.grid(
    refl = seq_along(reflec),
    gamma = gamma,
    omega = if (tno) omega else omega[1L],
    KEEP.OUT.ATTRS   = FALSE,
    stringsAsFactors = FALSE
  )

  # extract residuals metrics from Tb estimates
  get_metrics <- function(refl_idx, g, o) {
    res <- retvod::estTb(
      tbH   = tbH,
      tbV   = tbV,
      fH    = reflec[[refl_idx]]$fH,
      fV    = reflec[[refl_idx]]$fV,
      gamma = g,
      Tair  = Tair,
      Tsoil = Tsoil,
      omega = o
    )
   out <- c(
      pred_tbH = res$pred_tbH,
      pred_tbV = res$pred_tbV,
      cf_total = res$residuals$totaltb,
      cf_tbH   = res$residuals$tbH,
      cf_tbV   = res$residuals$tbV
    )
    return(out)
  }

  # 3) Compute metrics for every row of param_grid
  metrics_mat <- t( mapply(
    FUN      = get_metrics,
    refl_idx = grid$refl,
    g        = grid$gamma,
    o        = grid$omega,
    SIMPLIFY = TRUE
  ) )

  # 4) Combine into one results data.frame
  results_df <- cbind(grid, as.data.frame(metrics_mat))

  # 5) Identify the best (lowest cf_total)
  best_row <- results_df[ which.min(results_df$cf_total), ]

  # 6) Package outputs to mirror your original list structure
  out <- list(
    best = best_row,
    min_cf_index = c(best_row$refl, best_row$gamma, best_row$omega),
    cf_tb        = best_row$cf_total,
    pred_tbH     = best_row$pred_tbH,
    pred_tbV     = best_row$pred_tbV,
    cf_tbH       = best_row$cf_tbH,
    cf_tbV       = best_row$cf_tbV,
    reflec_best  = reflec[[best_row$refl]],
    gamma_best   = best_row$gamma,
    omega_best   = best_row$omega
  )

  # attach full cost grid if requested
  if (mat) {
    out$cf_mat <- results_df
  }

  return(out)
}


#' # Helper for the getting residuals
#'
#' #' Get residuals in grid across a range of gamma and omega.
#' #'
#' #' @param reflec reflectivities
#' #' @param g gamma range
#' #' @param o omega range
#' #'
#' #' @returns A matrix of residuals for each gamma and omega combination.
#' #' @export
#'
#' get_residuals <- function(reflec, g, o) {
#'
#'   out <- retvod::estTb(
#'     tbH = tbH,
#'     tbV = tbV,
#'     fH = reflec[[1]]$fH,
#'     fV = reflec[[1]]$fV,
#'     gamma = g,
#'     Tair = Tair,
#'     Tsoil = Tsoil,
#'     omega = o
#'   )$residuals$totaltb
#'
#'   return(out)
#' }
