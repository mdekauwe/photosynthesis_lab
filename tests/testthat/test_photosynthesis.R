library(testthat)

test_that("arrh returns expected value at 298K", {
  expect_equal(arrh(100, 50000, 298.15), 100)
})


test_that("peaked_arrh returns lower value at high temperature", {
  val1 <- peaked_arrh(100, 50000, 298.15, 650, 200000)
  val2 <- peaked_arrh(100, 50000, 313.15, 650, 200000)
  expect_lt(val2, val1)
})

test_that("calc_michaelis_menten_constants returns numeric", {
  p <- list(Kc25 = 404.9, Ec = 79430, Ko25 = 278.4, Eo = 36380, Oi = 210)
  val <- calc_michaelis_menten_constants(p, 298.15)
  expect_true(is.numeric(val))
})

test_that("calc_stomatal_coeff returns finite positive", {
  p <- list(g1 = 4, g0 = 0)
  val <- calc_stomatal_coeff(p, 400, 1)
  expect_true(all(val > 0))
})

test_that("calc_photosynthesis returns expected structure", {
  p <- list(
    Kc25 = 404.9, Ec = 79430, Ko25 = 278.4, Eo = 36380, Oi = 210,
    gamstar25 = 42.75, Eag = 37830,
    Vcmax25 = 40, Eav = 51560, deltaSv = 629.26, Hdv = 200000,
    Jmax25 = 66.8, Eaj = 43790, deltaSj = 644.4338, Hdj = 200000,
    theta_J = 0.7, alpha = 0.24, g0 = 1e-9, g1 = 4
  )
  result <- calc_photosynthesis(p, 298.15, 1500, 400, 1.5)
  expect_true(is.list(result))
  expect_named(result, c("An", "Ac", "Aj", "gsc", "Vcmax", "Cic", "Rd"))
  expect_true(all(sapply(result, is.numeric)))
})

