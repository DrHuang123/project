#' Perform Analysis on Gene Expression Data
#'
#' @description This function performs normalization and hypothesis testing on rhythmic gene expression data.
#'
#' @param data A numeric matrix of gene expression counts (rows: genes, columns: samples).
#' @param time A numeric vector indicating time points for each sample.
#' @param ncond An integer specifying the number of experimental conditions.
#' @param period A numeric value indicating the periodicity (default is 24).
#' @param verbose A logical value indicating whether to print progress messages (default is TRUE).
#' @return A list with the following components:
#'   \describe{
#'     \item{TestStatistics}{A numeric vector of test statistics for each gene.}
#'     \item{PValues}{A numeric vector of p-values for each gene.}
#'     \item{AdjustedPValues}{A numeric vector of adjusted p-values (Benjamini-Hochberg).}
#'   }
#' @examples
#' data(ExampleData)
#' data <- ExampleData[["CountData"]]
#' group <- ExampleData[["group"]]
#' time <- ExampleData[["time"]]
#' result <- BayesianRhythmicAnalysis(data, time, group, ncond = 2, period = 24, verbose = TRUE)
#' @export
BayesianRhythmicAnalysis <- function(data, time, group, ncond, period = 24, verbose = TRUE) {

  if (verbose) message("Checking DESeq2 package...")
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("The DESeq2 package is required but not installed.")
  }
  library(DESeq2)

  if (verbose) message("Filtering count data...")
  data <- as.matrix(data)
  countData <- data[rowSums(data) != 0, ]

  if (verbose) message("Creating DESeqDataSet...")
  group <- as.factor(group)
  colData <- data.frame(row.names = colnames(countData), group= group, time= time)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countData,
    colData   = colData,
    design    = ~ group
  )

  if (verbose) message("Calculating size factors...")
  dds <- DESeq2::estimateSizeFactors(dds)

  if (verbose) message("Normalizing data...")
  normalized_counts <- counts(dds, normalized = TRUE)
  mydata <- normalized_counts

  if (verbose) message("Building design matrix...")
  mymat <- matrix(
    c(sapply(time, function(t) c(1, cos(2*pi*t/period), sin(2*pi*t/period))),
      recursive = TRUE),
    ncol = 3, byrow = TRUE
  )

  if (verbose) message("Constructing hypothesis test matrix...")
  ncmymat <- ncol(mymat)
  nrmymat <- nrow(mymat)

  Lmat <- matrix(0, nrow = (ncond - 1) * 2, ncol = ncond * ncmymat)
  row_index <- 1
  for (ix in 2:ncond) {

    Lmat[row_index, 2] <- 1
    Lmat[row_index, (ix - 1)*ncmymat + 2] <- -1
    row_index <- row_index + 1

    Lmat[row_index, 3] <- 1
    Lmat[row_index, (ix - 1)*ncmymat + 3] <- -1
    row_index <- row_index + 1
  }


  if (verbose) message("Preparing variables and estimating error variance...")
  p <- nrow(mydata)
  store_var <- matrix(0, ncol = ncond, nrow = p)
  caps <- length(time)/(ncond*length(unique(time)))
  capi <- length(unique(time))

  xtx <- t(mymat) %*% mymat
  xtxinv <- solve(xtx)
  xtxinvxt <- xtxinv %*% t(mymat)

  nnull = 1000
  alpha = 0.05


  gamma_hat = betavec = nullgamma_hat = matrix(0, nrow=ncond, ncol=ncmymat)
  y = matrix(0, nrow=ncond, ncol=nrmymat)
  var_error = nullvar_error = rep(0, ncond)
  varcov = nullvarcov = matrix(0, nrow=(ncmymat*ncond), ncol=(ncmymat*ncond))

  yList <- vector("list", p)


  for (i in 1:p) {
    for (i1 in 1:ncond) {
      i2 <- (i1 - 1)*caps*capi + 1
      i3 <- i1*caps*capi
      y[i1, ] <- mydata[i, i2:i3]

      gamma_hat[i1, ] <- xtxinvxt %*% y[i1, ]
      residuals <- y[i1, ] - mymat %*% gamma_hat[i1, ]
      var_error[i1] <- t(residuals) %*% residuals / (nrmymat - 3)
    }
    store_var[i, ] <- var_error
    yList[[i]] <- y
  }

  if (verbose) message("Estimating prior distribution shape and scale parameters of the error variance...")
  sigma2gc <- matrix(0, nrow=p, ncol=ncond)
  psic0 <- rep(0, ncond)
  psic1 <- rep(0, ncond)

  for (i1 in 1:ncond) {
    sigma2gc[, i1] <- store_var[, i1]
    psic0[i1] <- (2 + (mean(sigma2gc[, i1]^2)/(mean(sigma2gc[, i1])^2) - 1)^(-1))
    psic1[i1] <- 1 / ((psic0[i1]-1)*mean(sigma2gc[, i1]))
  }


  if (verbose) message("Fitting the model and computing null distribution...")

  myfunc1 <- function(i, yMatrix) {

    gamma_hat <- betavec <- nullgamma_hat <- matrix(0, nrow=ncond, ncol=ncmymat)
    nully     <- matrix(0, nrow=ncond, ncol=nrmymat)
    varcov    <- nullvarcov <- matrix(0, nrow=(ncmymat*ncond), ncol=(ncmymat*ncond))

    for (i1 in 1:ncond) {
      gamma_hat[i1, ] <- xtxinvxt %*% yMatrix[i1, ]
      q <- yMatrix[i1, ] - mymat %*% gamma_hat[i1, ]

      term1 <- as.numeric(1 + 0.5*psic1[i1]*sum(q^2))
      term11 <- caps*capi/2 + psic0[i1]
      term12 <- term11 * (
        -psic1[i1]*xtx/term1 +
          psic1[i1]^2 * t(mymat)%*%q%*%t(q)%*%mymat/(term1^2)
      )
      idx_start <- (i1-1)*ncmymat + 1
      idx_end   <- i1*ncmymat
      varcov[idx_start:idx_end, idx_start:idx_end] <- as.matrix(-solve(term12))
    }

    theta_hat <- as.vector(t(gamma_hat))
    Lcov <- solve(Lmat %*% varcov %*% t(Lmat))
    tssec2 <- as.numeric(t(Lmat %*% theta_hat) %*% Lcov %*% (Lmat %*% theta_hat))

    # estimation of the null distribution
    betavec[, 1] <- gamma_hat[, 1]
    tempo1 <- apply(gamma_hat[, -1], 2, mean)
    for (i1 in 1:ncond) {
      betavec[i1, -1] <- tempo1
    }

    myfunull <- function(itn) {
      set.seed(i*itn)
      for (i1 in 1:ncond) {
        nully[i1, ] <- mymat %*% betavec[i1, ] +
          rnorm(caps*capi)*sqrt(store_var[i, i1])

        nullgamma_hat[i1, ] <- xtxinvxt %*% nully[i1, ]
        nullq <- nully[i1, ] - mymat %*% nullgamma_hat[i1, ]
        nullqtq <- sum(nullq^2)

        nullterm1 <- as.numeric(1 + 0.5*psic1[i1]*nullqtq)
        nullterm11 <- caps*capi/2 + psic0[i1]
        nullterm12 <- nullterm11*(
          -psic1[i1]*xtx/nullterm1 +
            psic1[i1]^2 * t(mymat)%*%q%*%t(q)%*%mymat/(nullterm1^2)
        )
        idx_start <- (i1-1)*ncmymat + 1
        idx_end   <- i1*ncmymat
        nullvarcov[idx_start:idx_end, idx_start:idx_end] <- as.matrix(-solve(nullterm12))
      }
      nulltheta_hat <- as.vector(t(nullgamma_hat))
      nullLcov <- solve(Lmat %*% nullvarcov %*% t(Lmat))
      tempo2 <- Lmat %*% nulltheta_hat
      nullresults <- as.numeric(t(tempo2) %*% nullLcov %*% tempo2)
      return(nullresults)
    }
    nullresults <- as.numeric(lapply(1:nnull, myfunull))
    pvaluesec2  <- length(which(nullresults > tssec2)) / nnull

    df_this_gene <- data.frame(
      Gene      = i,
      Condition = 1:ncond,
      Alpha     = gamma_hat[,1],
      BetaCos   = gamma_hat[,2],
      BetaSin   = gamma_hat[,3]
    )
    attr(df_this_gene, "TestStat") <- tssec2
    attr(df_this_gene, "PValue")   <- pvaluesec2
    return(df_this_gene)
  }


  if (verbose) message("Running parallel analysis...")
  library(doParallel)
  registerDoParallel(cores = parallel::detectCores() - 1)

  out_list <- foreach(it = 1:p) %dopar% {
    myfunc1(it, yList[[it]])
  }

  if (verbose) message("Gathering test statistics and p-values...")
  teststat   <- sapply(out_list, function(x) attr(x, "TestStat"))
  pvaluesec2 <- sapply(out_list, function(x) attr(x, "PValue"))

  res_df <- do.call(rbind, out_list)

  if (verbose) message("Adjusting p-values...")
  adj_pvalues <- p.adjust(pvaluesec2, method = "BH")

  if (verbose) message("Adding gene names...")
  gene_names <- rownames(data)
  res_df$GeneName <- rep(gene_names[1:p], each = ncond)

  if (verbose) message("Returning final results...")
  return(list(
    TestStatistics = teststat,
    PValues        = pvaluesec2,
    AdjustedPValues= adj_pvalues,
    FittedParams   = res_df
  ))
}
