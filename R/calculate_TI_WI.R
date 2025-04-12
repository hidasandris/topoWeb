#' Calculate TI/WI index
#'
#' This function calculates the TI and WI index for networks. Input should be an igraph object or
#' a data frame without header, the columns should be ordered: node name of species A,
#' node name of species B, weight (if exists).
#' Symmetrization_method is optional, can be either 'sum' (default), 'average', 'max', 'min' or
#' 'difference.' Multilink edges will be summed.
#'
#' @param input_data Input igraph object or data frame
#' @param steps TI steps to compute
#' @param symmetrization_method How to combine TI values
#' @param asymmetry Whether to return an input for asymmetry calculation
#' @returns A data frame with node ids and TI, WI (if there are weights) index.
#' @export
calculate_TI_WI <- function(input_data, steps, symmetrization_method = "sum", asymmetry = F) {
  numstep <- steps
  symtype <- symmetrization_method

  if (class(input_data) == "igraph") {
    if (!require(igraph)) stop("igraph package not installed!")
    input_data_edgelist <- igraph::as_edgelist(input_data)

    if (is.null(E(input_data)$weight)) {
      input_data_edgelist <- data.frame(
        node1 = input_data_edgelist[, 1],
        node2 = input_data_edgelist[, 2]
      )
      input_data <- input_data_edgelist[c("node1", "node2")]
    } else {
      input_data_edgelist <- data.frame(
        node1 = input_data_edgelist[, 1],
        node2 = input_data_edgelist[, 2],
        weight = E(input_data)$weight
      )
      input_data <- input_data_edgelist[c("node1", "node2", "weight")]
    }
  }

  if (class(input_data) != "data.frame") input_data <- as.data.frame(input_data)

  if (!dim(input_data)[2] %in% c(2, 3)) stop("Wrong input format")

  if (dim(input_data)[2] == 2) {
    if (is_character(input_data[1, 1])) {
      names(input_data) <- c("V1", "V2")
    }
  } else {
    if (is_character(input_data[1, 1])) {
      names(input_data) <- c("V1", "V2", "V3")
    }

    if (any(duplicated(input_data[, 1:2]))) {
      input_data <- aggregate(V3 ~ V1 + V2, data = input_data, FUN = sum)
      warning("Multilink edges were summed.")
    }
  }

  node_id <- unique(c(input_data[, 1], input_data[, 2]))
  numnode <- length(node_id)
  mx <- matrix(rep(0, numnode^2), nrow = numnode, ncol = numnode)
  rownames(mx) <- node_id
  colnames(mx) <- node_id
  for (i in 1:length(input_data[, 1]))
  {
    mx[as.character(input_data[i, 1]), as.character(input_data[i, 2])] <- 1
    mx[as.character(input_data[i, 2]), as.character(input_data[i, 1])] <- 1
  }
  TI <- mx
  for (i in 1:numnode)
  {
    for (j in 1:numnode) {
      TI[j, i] <- mx[j, i] / sum(mx[, i])
    }
  }
  SI <- matrix(rep(0, numnode^2), nrow = numnode, ncol = numnode)
  CI <- diag(numnode)
  for (i in 1:numstep)
  {
    CI <- CI %*% TI
    SI <- SI + CI
  }
  TI_index <- numeric(numnode)
  for (i in 1:numnode) {
    TI_index[i] <- sum(SI[i, ]) / numstep / numnode
  }
  resu <- data.frame(node_id, TI_index)


  if (dim(input_data)[2] == 2) {
    warning("No weights, computing only TI index.")
    if (asymmetry) {
      return(list(SI = SI, node_id = node_id))
    }
    return(resu)
  }

  if (!is.numeric(input_data[[3]])) stop("Weights are not of numeric type")

  if (any(input_data[, 3] < 0)) {
    input_data[, 3] <- abs(input_data[, 3])
    warning("Negative values detected. Absolute values taken.")
  }

  mxw <- matrix(rep(0, numnode^2), nrow = numnode, ncol = numnode)
  rownames(mxw) <- node_id
  colnames(mxw) <- node_id
  for (i in 1:length(input_data[, 1]))
  {
    mxw[as.character(input_data[i, 1]), as.character(input_data[i, 2])] <-
      input_data[i, 3]
  }
  for (i in 1:numnode) {
    for (j in 1:numnode) {
      if (i < j) {
        if ((mxw[i, j] == 0) ||
          (mxw[j, i] == 0)) {
          va <- mxw[i, j] + mxw[j, i]
          mxw[i, j] <- va
          mxw[j, i] <- va
        } else {
          if (symtype == "sum") {
            va <- mxw[i, j] + mxw[j, i]
          }
          if (symtype == "average") {
            va <- (mxw[i, j] + mxw[j, i]) / 2
          }
          if (symtype == "max") {
            va <- max(mxw[i, j], mxw[j, i])
          }
          if (symtype == "min") {
            va <- min(mxw[i, j], mxw[j, i])
          }
          if (symtype == "difference") {
            va <- abs(mxw[i, j] - mxw[j, i])
          }
          mxw[i, j] <- va
          mxw[j, i] <- va
        }
      }
    }
  }
  TI <- mxw
  for (i in 1:numnode)
  {
    for (j in 1:numnode) {
      TI[j, i] <- mxw[j, i] / sum(mxw[, i])
    }
  }
  SI <- matrix(rep(0, numnode^2), nrow = numnode, ncol = numnode)
  CI <- diag(numnode)
  for (i in 1:numstep)
  {
    CI <- CI %*% TI
    SI <- SI + CI
  }
  WI_index <- numeric(numnode)
  for (i in 1:numnode) {
    WI_index[i] <- sum(SI[i, ]) / numstep / numnode
  }
  resuw <- data.frame(node_id, WI_index)

  result <- cbind(resu, WI_index)

  if (asymmetry) {
    return(list(SI = SI, node_id = node_id))
  }
  return(result)
}
