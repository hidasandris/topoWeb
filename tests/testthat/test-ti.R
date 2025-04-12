toy_data_TI3 <- toy_data |>
  calculate_TI_WI(3)

toy_data_TI3_igraph <- toy_data |>
  igraph::graph_from_data_frame() |>
  calculate_TI_WI(3)

toy_data_TI3_expected <- data.frame(
  node_id = LETTERS[1:8],
  TI_index = c(
    0.05227623,
    0.18798225,
    0.09958526,
    0.15186150,
    0.23389275,
    0.05227623,
    0.17481674,
    0.04730903
  )
)

test_that("TI3 from data.frame ok", {
  expect_equal(toy_data_TI3, toy_data_TI3_expected, tolerance = 0.0001)
})

test_that("TI3 from igraph ok", {
  expect_equal(toy_data_TI3_igraph, toy_data_TI3_expected, tolerance = 0.0001)
})
