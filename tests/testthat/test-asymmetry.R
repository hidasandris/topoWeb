toy_data_asym <- toy_data |>
  calculate_asymmetry(filter = F)

test_that("Number of asymmetry values ok", {
  expect_equal(nrow(toy_data_asym), 24)
})


toy_data_asym_default_filter <- toy_data |>
  calculate_asymmetry()

test_that("Number of default filtered asymmetry values ok", {
  expect_equal(nrow(toy_data_asym_default_filter), 0)
})


toy_data_asym_filter_0.2 <- toy_data |>
  calculate_asymmetry(threshold = 20)

toy_data_asym_filter_0.2_expected <- data.frame(
  from = c("B", "E", "G"),
  to = c("A", "H", "F"),
  weight = c(44.75309, 48.17708, 42.43827)
)

test_that("Asymmetry values ok", {
  expect_equal(toy_data_asym_filter_0.2,
               toy_data_asym_filter_0.2_expected,
               tolerance = 0.0001)
})
