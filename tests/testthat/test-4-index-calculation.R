test_that("get_index(), get_index_sims(), and get_cog() work", {
  local_edition(3)
  skip_on_ci()
  skip_on_cran()
  skip_if_not_installed("INLA")
  pcod_spde <- make_mesh(pcod, c("X", "Y"), cutoff = 20)
  m <- sdmTMB(
    data = pcod,
    formula = density ~ 0 + as.factor(year),
    time = "year", mesh = pcod_spde, family = tweedie(link = "log")
  )
  # expect_snapshot(m)
  predictions <- predict(m, newdata = qcs_grid, return_tmb_object = TRUE)

  p <- predict(m, newdata = qcs_grid, return_tmb_object = FALSE)
  expect_error(get_index(p), regexp = "return_tmb_object")

  ind <- get_index(predictions, bias_correct = FALSE)
  expect_equal(class(ind), "data.frame")

  ind_corrected <- get_index(predictions, bias_correct = TRUE)
  cached <- c(263131.7738, 390965.6367, 432988.7341, 119184.5699, 192717.7117,
    338393.0779, 341000.346, 405631.4672, 195346.191)
  expect_lt(sum(abs(ind_corrected$est - cached) / cached), 1e-03)
  expect_gt(mean(ind_corrected$est - ind$est), 0)

  cached <- c(209276.6742, 316827.7388, 355482.8274, 91541.2387, 150187.3753,
    274230.3944, 274541.7833, 326960.9376, 151502.4236)
  expect_lt(sum(abs(ind$est - cached) / cached), 1e-03)

  cached <- c(151512.8161, 244649.3652, 274277.5841, 66923.457, 109088.6819,
    212808.0659, 204421.8674, 241466.6371, 110227.6953)
  expect_lt(sum(abs(ind$lwr - cached) / cached), 1e-03)

  set.seed(1)
  pred_sim <- predict(m, nsim = 2000L)
  ind_sim <- get_index_sims(pred_sim)
  expect_gt(cor(ind_sim$est, ind$est), 0.9)
  expect_gt(cor(ind_sim$lwr, ind$lwr), 0.9)
  expect_gt(cor(ind_sim$upr, ind$upr), 0.9)
  # sims mimics bias corrected index, which would be higher:
  expect_gt(mean(ind$est), mean(ind_sim$est))

  cog <- get_cog(predictions, bias_correct = FALSE)
  expect_equal(class(cog), "data.frame")
  cached <- c(462.0865, 481.3353, 471.6714, 481.9013, 485.464, 469.7996,
    475.9824, 457.4335, 463.3973, 5756.9909, 5729.2348, 5760.9857,
    5735.0945, 5726.1765, 5745.2497, 5744.2829, 5757.1354, 5755.7778
  )
  expect_lt(sum(abs(cog$est - cached) / cached), 1e-04)

  # missing time:
  .qcs_grid <- subset(qcs_grid, year != 2015)
  expect_error(p <- predict(m, newdata = .qcs_grid, return_tmb_object = TRUE), regexp = "time")
})

