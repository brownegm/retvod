test_that("Check errors", {
  # check that the lengths are the same.
  expect_error(
    retVOD(
      tbH = seq(100, 200, length.out = 10),
      tbV = seq(100, 200, length.out = 10),
      smc = seq(0.5, 0.9, length.out = 9), vod = 1, Tair = 1, Tsoil = 1,
      omega = 1, roughness = 1, inc_angle = 1
    ),
    "tbH, tbV, and smc lengths differ."
  )
})

test_that("Output matches manual single measurement solution", {
  v <- c(
    268.4722473, 268.003972, 267.1787049,
    266.446253, 266.4795018, 266.1335705,
    265.6720042
  )
  h <- c(
    258.957763, 258.4377463, 257.8380787,
    257.2369798, 257.1168902, 256.7121893,
    256.5970527
  )

  air <- c(
    291.32, 291.3, 290.2,
    289.24, 288.76, 288.12,
    287.88
  )

  soil <- air * 0.90

  set.seed(2)
  sm <- retvod:::gen_sin(1000, rangeL = 0.2, rangeH = 0.45) |> sample(size = 7)
  vod <- seq(0,2, length.out = 7)
  inc_angle <- 40
  ## calculate gamma for each VOD test value
  gamma <- exp(-vod / cos(inc_angle * (pi / 180)))
  clay_frac=0.232
  ## calculate epsilon (dielectric) and reflectivitys for each value of soil moisture
  eps_list <- sapply(sm, \(s) mironov(1.4e9, s, clay_frac)$dielectric)
  reflecs <- sapply(eps_list, \(e) fresnelr(eps = e, theta = inc_angle, h = 0.1), simplify = F)

  sol2 <- solveSmVod(
    reflecs[4], gamma = gamma, tbH = h[4], tbV = v[4],
    Tair = air[4], Tsoil = soil[4],
    omega = 0.05, #clay_frac = 0.232,
    mat = T
  )

  retrieval <- retVOD(h,v,
                      smc = sm, vod = vod,
                      Tair = air, Tsoil = soil, cf = 0.232,
                      omega = 0.05, roughness = 0.1, inc_angle = 40)

  expect_identical(sm, retrieval$smEst) # soil moisture input and output are same
  H <- retrieval[4,"rH"] == sol2$reflec_best$fH
  V <- retrieval[4,"rV"] == sol2$reflec_best$fV
  expect_true(H && V) # reflectivity values are same as in solution
  #expect_equal(retrieval[4,c("rH","rV")], sol2$reflec_best) # output from solution is same as retrieved values
  expect_equal(retrieval$cfEst[4],sol2$cf_tb)
  expect_equal(retrieval$tbVpred[4], sol2$pred_tbV)
  expect_equal(retrieval$tbHpred[4], sol2$pred_tbH)
  expect_equal(retrieval$tbHcost[4], sol2$cf_tbH)
  expect_equal(retrieval$tbVcost[4], sol2$cf_tbV)
  expect_s3_class(retrieval,"retVOD")

#
#   vodestimates <- retrieval$vodEst
#
#   retrieval2 <- retVOD(h,v,
#                       smc = seq(0.2, 0.45, length.out = 7), vod = vodestimates,
#                       Tair = air, Tsoil = soil, cf = 0.232,
#                       omega= 0.05, roughness = 0.1, inc_angle = 40)

})


test_that("DUAL Omega and Tau: Output matches manual single measurement solution", {
  v <- c(
    268.4722473, 268.003972, 267.1787049,
    266.446253, 266.4795018, 266.1335705,
    265.6720042
  )
  h <- c(
    258.957763, 258.4377463, 257.8380787,
    257.2369798, 257.1168902, 256.7121893,
    256.5970527
  )

  air <- c(
    291.32, 291.3, 290.2,
    289.24, 288.76, 288.12,
    287.88
  )

  soil <- air * 0.90

  set.seed(2)
  sm <- retvod:::gen_sin(1000, rangeL = 0.2, rangeH = 0.45) |> sample(size = 7)
  vod <- seq(0.2,3, by = 0.1)
  inc_angle <- 40
  ## calculate gamma for each VOD test value
  gamma <- exp(-vod / cos(inc_angle * (pi / 180)))
  clay_frac = 0.232
  ## calculate epsilon (dielectric) and reflectivitys for each value of soil moisture
  eps_list <- sapply(sm, \(s) mironov(1.4e9, s, clay_frac)$dielectric)
  reflecs <- sapply(eps_list, \(e) fresnelr(eps = e, theta = inc_angle, h = 0.1), simplify = F)
  # log_omegas <- seq(log(0.01), log(1), length.out = 200)
  # omega_range <- exp(log_omegas)

  omega_range <- seq(0.1, 0.5, by = 0.01)

  sol2 <- solveSmVod(
    reflecs[4], gamma = gamma, tbH = h[4], tbV = v[4],
    Tair = air[4], Tsoil = soil[4],
    omega = omega_range, tno = T,
    mat = T)

  plot(sol2$cf_mat$cf_total~sol2$cf_mat$omega)

  retrieval <- retVOD(h,v,
                      smc = sm, vod = vod,
                      Tair = air, Tsoil = soil, cf = 0.232,
                      omega = omega_range, roughness = 0.1, inc_angle = 40, tno = T)

  expect_identical(sm, retrieval$smEst) # soil moisture input and output are same
  H <- retrieval[4,"rH"] == sol2$reflec_best$fH
  V <- retrieval[4,"rV"] == sol2$reflec_best$fV
  expect_true(H && V) # reflectivity values are same as in solution
  #expect_equal(retrieval[4,c("rH","rV")], sol2$reflec_best) # output from solution is same as retrieved values
  expect_equal(retrieval$cfEst[4],sol2$cf_tb)
  expect_equal(retrieval$tbVpred[4], sol2$pred_tbV)
  expect_equal(retrieval$tbHpred[4], sol2$pred_tbH)
  expect_equal(retrieval$tbHcost[4], sol2$cf_tbH)
  expect_equal(retrieval$tbVcost[4], sol2$cf_tbV)
  expect_s3_class(retrieval,"retVOD")

})
