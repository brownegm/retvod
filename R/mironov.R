#' mironov
#'
#' Calculate the dielectric constant of soil
#' @name mironov
#' @param f Frequency
#' @param smv Soil moisture volume
#' @param cf Clay fraction of soil
#' @return A list of the real and imaginary parts of the dielectric constant
#' @export
#'
mironov <- function(f, smv, cf) {

  fcheck <- f <= 1e9
  #smcheck <- smv <= 0 | smv >= 1
  cfcheck <- cf <= 0 | cf >= 1

  if (fcheck | cfcheck) {# smcheck |
    cat("\nfcheck =", fcheck,"\n",
              #"smcheck =", smcheck, "\n",
              "cfcheck =", cfcheck, "\n")
    stop("mironov: Check input---Data falls outside of acceptable range.")
  }
  # Constants
  eps_0 <- 8.854 * 10^(-12) # static dielectric constant of water
  eps_winf <- 4.9 # high frequency limit of dielectric constant of water
  fHz <- f # frequency in hertz

  # Initializing the GRMDM spectroscopic parameters with clay fraction
  # RI & NAC of dry soils
  znd <- 1.634 - 0.539 * cf + 0.2748 * cf^2
  zkd <- 0.03952 - 0.04038 * cf

  # Maximum bound water fraction
  zxmvt <- 0.02863 + 0.30673 * cf
  # zxmvt = 0.02863 + 0.0030673 * cf

  # Bound water parameters
  zep0b <- 79.8 - 85.4 * cf + 32.7 * cf^2
  ztaub <- 1.062 * 10^(-11) + 3.450 * 10^(-12) * cf
  zsigmab <- 0.3112 + 0.467 * cf

  # Unbound (free) water parameters
  zep0u <- 100
  ztauu <- 8.5 * 10^(-12)
  zsigmau <- 0.3631 + 1.217 * cf

  # Computation of epsilon water (bound & unbound)
  zcxb <- (zep0b - eps_winf) / (1 + (2 * pi * fHz * ztaub)^2)
  zepwbx <- eps_winf + zcxb
  zepwby <- zcxb * (2 * pi * fHz * ztaub) + zsigmab / (2 * pi * eps_0 * fHz)
  zcxu <- (zep0u - eps_winf) / (1 + (2 * pi * fHz * ztauu)^2)
  zepwux <- eps_winf + zcxu
  zepwuy <- zcxu * (2 * pi * fHz * ztauu) + zsigmau / (2 * pi * eps_0 * fHz)

  # Computation of refractive index of water (bound & unbound)
  znb <- sqrt(sqrt(zepwbx^2 + zepwby^2) + zepwbx) / sqrt(2)
  zkb <- sqrt(sqrt(zepwbx^2 + zepwby^2) - zepwbx) / sqrt(2)
  znu <- sqrt(sqrt(zepwux^2 + zepwuy^2) + zepwux) / sqrt(2)
  zku <- sqrt(sqrt(zepwux^2 + zepwuy^2) - zepwux) / sqrt(2)

  # Computation of soil refractive index (nm & km): xmv can be a vector
  zxmvt2 <- min(smv, zxmvt)
  zflag <- smv >= zxmvt
  znm <- znd + (znb - 1) * zxmvt2 + (znu - 1) * (smv - zxmvt) * zflag
  zkm <- zkd + zkb * zxmvt2 + zku * (smv - zxmvt) * zflag

  # Computation of soil dielectric constant:
  zepmx <- znm^2 - zkm^2
  zepmy <- znm * zkm * 2
  # eps = zepmx + i*zepmy
  eps <- complex(real = zepmx, imaginary = -1 * zepmy)

  # return the real and imaginary parts of the dielectric constant
  return(structure(list(zepmx, zepmy, eps),
    names = c("real", "imaginary", "dielectric")
  ))
}
