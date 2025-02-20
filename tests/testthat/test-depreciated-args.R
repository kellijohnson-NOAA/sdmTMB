test_that("Depreciated args work/throw warnings/stops", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")

  # fields
  expect_warning(m1 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    fields = "AR1",
    time = "year",
    mesh = pcod_mesh_2011,
    do_fit = FALSE
  ), regexp = "fields")
  m2 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    spatiotemporal = "AR1",
    time = "year",
    mesh = pcod_mesh_2011,
    do_fit = FALSE
  )
  expect_warning(m3 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    ar1_fields = TRUE,
    time = "year",
    mesh = pcod_mesh_2011,
    do_fit = FALSE
  ), "ar1_fields")
  expect_identical(m1$tmb_data, m2$tmb_data)
  expect_identical(m1$tmb_data, m3$tmb_data)

  # spatial_only
  expect_warning(m1 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    spatial_only = TRUE,
    time = "year",
    mesh = pcod_mesh_2011,
    do_fit = FALSE
  ), regexp = "spatial_only")
  m2 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    spatiotemporal = "off",
    time = "year",
    mesh = pcod_mesh_2011,
    do_fit = FALSE
  )
  expect_identical(m1$tmb_data, m2$tmb_data)
  m3 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    time = NULL,
    mesh = pcod_mesh_2011,
    do_fit = FALSE
  )
  expect_identical(m1$tmb_data$spatial_only, m3$tmb_data$spatial_only)
  expect_identical(m3$tmb_data$n_t, 1L)
  expect_identical(m1$tmb_data$n_t, 4L)

  # include_spatial
  expect_warning(m1 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    include_spatial = FALSE,
    time = "year",
    mesh = pcod_mesh_2011,
    do_fit = FALSE
  ), regexp = "include_spatial")
  m2 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    spatial = "off",
    time = "year",
    mesh = pcod_mesh_2011,
    do_fit = FALSE
  )
  expect_identical(m1$tmb_data, m2$tmb_data)
  m2 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    spatial = FALSE,
    time = "year",
    mesh = pcod_mesh_2011,
    do_fit = FALSE
  )
  expect_identical(m1$tmb_data, m2$tmb_data)
})

test_that("spatial = off and spatiotemporal = off triggers map_rf", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  m1 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    spatial = "off",
    spatiotemporal = "off",
    mesh = pcod_mesh_2011
  )
  m2 <- sdmTMB(
    formula = density ~ 1,
    data = pcod_2011,
    control = sdmTMBcontrol(map_rf = TRUE),
    mesh = pcod_mesh_2011
  )
  expect_equal(m1$model$par, m2$model$par)
})

test_that("spde/mesh args work", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  expect_warning(
    m1 <- sdmTMB(density ~ 1, data = pcod_2011, spde = pcod_mesh_2011, do_fit = FALSE),
    "mesh"
  )
  m2 <- sdmTMB(density ~ 1, data = pcod_2011, mesh = pcod_mesh_2011, do_fit = FALSE)
  expect_identical(m1$tmb_data, m2$tmb_data)
})
