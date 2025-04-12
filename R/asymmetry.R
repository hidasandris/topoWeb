#' Calculate asymmetry values
#'
#' This function calculates interaction asymmetry based on TI index for networks.
#' Input should be an igraph object or a data frame without header,
#' the columns should be ordered: node name of species A,
#' node name of species B, weight (if exists).
#' If filter == TRUE (default) then top 1% of asymmetric links will be returned. Other absolute
#' threshold values can be supplied.
#'
#' @param input_data Input igraph object or data frame
#' @param threshold Threshold to use, default is top 1% of asymmetric links
#' @param TI_steps TI steps to compute
#' @param filter Filter the most asymmetric links
#' @returns A data frame with links indicating direction of asymmetric effect and it's value.
#' @importFrom rlang .data
#' @export
calculate_asymmetry <- function(input_data, threshold = NULL, TI_steps = 3, filter = TRUE) {
  TI_list <- suppressWarnings(calculate_TI_WI(input_data, TI_steps, asymmetry = TRUE))

  node_id <- TI_list$node_id
  numnode <- length(node_id)
  SI <- TI_list$SI / numnode / TI_steps * 1000

  AI <- SI - t(SI)
  row.names(AI) <- colnames(AI)

  AI_list_rows <- list()

  for (i in 1:numnode) {
    for (j in 1:numnode) {
      if (AI[i, j] <= 0) next
      AI_list_rows[[length(AI_list_rows) + 1]] <- data.frame(
        from = node_id[i], to = node_id[j], weight = AI[i, j]
      )
    }
  }

  AI_df <- dplyr::bind_rows(AI_list_rows)

  if (!filter) {
    return(AI_df)
  }

  if (is.null(threshold)) {
    asymmetry_filtered <- dplyr::slice_max(AI_df, .data$weight, prop = .01)
  } else {
    asymmetry_filtered <- AI_df[abs(AI_df$weight) >= threshold, ]
    row.names(asymmetry_filtered) <- NULL
  }

  return(asymmetry_filtered)
}
