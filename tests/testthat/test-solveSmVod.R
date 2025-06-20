test_that("Correct VOD values chosen", {
  # test values
  sm <- seq(0.2, 0.8, by = 0.2)
  vod_test <- seq(0.5, 1, by = 0.2)
  air <- 300
  soil <- 275
  omega <- 0.03
  inc_angle <- 40 # degrees
  clay_frac <- 0.232
  h <- 0.16
  tbH <- 280
  tbV <- 285
  sm <- sm[1]
  omega_range <- seq(0.01, 0.09, by = 0.01)
  ## calculate gamma for each VOD test value
  gamma <- exp(-vod_test / cos(inc_angle * (pi / 180)))

  ## calculate epsilon (dielectric) and reflectivitys for each value of soil moisture
  eps_list <- sapply(sm, \(s) mironov(1.4e9, s, clay_frac)$dielectric)
  reflecs <- lapply(eps_list, \(e) fresnelr(eps = e, theta = inc_angle, h))

  # run function
  sol1 <- solveSmVod(
    reflec = reflecs, gamma = gamma,
    tbH = tbH, tbV = tbV,
    Tair = air, Tsoil = soil,
    omega = omega_range,
    mat = T, tno = T
  )

  # check that the correct values are chosen
  expect_equal(reflecs[[sol1$min_cf_index[1]]], sol1$reflec_best)
  expect_true(sol1$gamma_best %in% gamma)
  expect_equal(min(sol1$cf_mat), sol1$cf_tb)

  # check specific indices
  expect_equal(reflecs[[1]], sol1$reflec_best)
  expect_equal(gamma[3], sol1$gamma_best)

  # check that the H and V cf sum is the same as what is reported for the total cf
  sum_cfHV <- sol1$cf_tbV + sol1$cf_tbH
  expect_equal(sum_cfHV, sol1$cf_tb)
})

test_that("Same values backward and forward", {
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
  sm <- retvod:::gen_sin(1000, rangeL = 0.2, rangeH = 0.35) |> sample(size = 7)
  vod <- sm * 2
  inc_angle <- 40
  clay_frac <- 0.3
  ## calculate gamma for each VOD test value
  gamma <- exp(-vod / cos(inc_angle * (pi / 180)))

  ## calculate epsilon (dielectric) and reflectivitys for each value of soil moisture
  eps_list <- sapply(sm, \(s) mironov(1.4e9, s, clay_frac)$dielectric)
  reflecs <- sapply(eps_list, \(e) fresnelr(eps = e, theta = inc_angle, h = 0.1), simplify = F)

  sol2forward <- solveSmVod(
    reflec = reflecs[4], gamma = gamma, tbH = h[4], tbV = v[4],
    Tair = air[4], Tsoil = soil[4],
    omega = 0.5,
    mat = T
  )

  expect_equal(min(sol2forward$cf_mat), sol2forward$cf_tb)
  expect_equal(sol2forward$reflec_best, reflecs[[4]])

  solved_gamma <- sol2forward$gamma_best
  input_sm <- sm[4]
  input_reflec <- reflecs[[4]]

  sol2_back <- solveSmVod(
    reflec = reflecs, gamma = solved_gamma, tbH = h[4], tbV = v[4],
    Tair = air[4], Tsoil = soil[4],
    omega = 0.5,
    mat = T
  )
  expect_identical(sol2forward[-c(1, 9)], sol2_back[-c(1, 9)])
})

