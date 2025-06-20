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
  num_r <- length(reflec) # number of dielectric values
  num_gamma <- length(gamma) # number of test VOD values
  num_omega <- length(omega) # number of omega values
  print(num_r)
  # array to hold residuals and predictions for each epsilon and gamma
  results_dim <- if (tno) {
    c(num_omega, num_gamma, 6)
  } else {
    c(num_r, num_gamma, 6)
}

  results <- array(NA,
    dim = results_dim, # rows, columns, calc values
    dimnames = list(seq(1,results_dim[1]),
                    seq(1,num_gamma),
                    c("pred_tbH", "pred_tbV", "cf_total", "cf_tbH", "cf_tbV", "omega"))
  )

  ## Compute cost function for all combinations of epsilon and VOD
  ## Q: What combo of smc and vod would the reflecs represent?
  ## For each eps(dielectric) and potential g(gamma,vod), calculate brightness temperatures
  ## and report the cost function (diff between calculated and obs Tbs)
  ## Doing a single variable retrieval right so one smc
  # if you want to do a retrieval with VOD and omega varying simulaneously
  if (tno) {
  for(e in seq_len(num_r)) { # loop through each reflec
    for (o in seq_len(num_omega)) {
      for (g in seq_len(num_gamma)) {
        omega <- omega[o] # get current omega value
        print(omega)
        result <- estTb(
          tbH = tbH,
          tbV = tbV,
          fH = reflec[[e]]$fH,
          fV = reflec[[e]]$fV,
          gamma = gamma[g],
          Tair = Tair,
          Tsoil = Tsoil,
          omega = omega
        )
        # store results
        results[o, g, ] <- c(
          result$pred_tbH,
          result$pred_tbV,
          result$residuals$totaltb,
          result$residuals$tbH,
          result$residuals$tbV,
          omega
        )
      }
    }
}
    } else { #traditional non varying omega

      for (e in seq_len(num_r)) {
        for (g in seq_len(num_gamma)) {
          result <- estTb(
            tbH = tbH,
            tbV = tbV,
            fH = reflec[[e]]$fH,
            fV = reflec[[e]]$fV,
            gamma = gamma[g],
            Tair = Tair,
            Tsoil = Tsoil,
            omega = omega
          )
          # store results
          results[e, g, ] <- c(
            result$pred_tbH,
            result$pred_tbV,
            result$residuals$totaltb,
            result$residuals$tbH,
            result$residuals$tbV
          )
        }
      }
    }

  min_index <- which(results[, , "cf_total"] == min(results[, , "cf_total"],
                                                    na.rm = TRUE), arr.ind = TRUE)
browser()
  if(num_r==1){
    best_row <- 1
    best_col <- min_index
  }else if (num_gamma==1){
    best_row <- min_index
    best_col <- 1
  }else{
    best_row <- min_index[1]
    best_col <- min_index[2]
  }

  output <- list(
    min_cf_index = c(best_row, best_col)|>unname(),
    cf_tb = results[best_row, best_col, "cf_total"]|>unname(),
    #epsilon = eps_list[best_row],
    pred_tbH = results[best_row, best_col, "pred_tbH"]|>unname(),
    pred_tbV = results[best_row, best_col, "pred_tbV"]|>unname(),
    cf_tbH = results[best_row, best_col, "cf_tbH"]|>unname(),
    cf_tbV = results[best_row, best_col, "cf_tbV"]|>unname(),
    reflec_best = reflec[[best_row]],
    gamma_best = gamma[best_col],
    cf_mat = if (mat == T) results[, , "cf_total"] else NA,
    omega = if (tno) results[best_row, best_col, "omega"] else NA
  )
  return(structure(output,
                   flag = if (length(best_row)>1) "Tie for lowest residuals found")
  )
}
