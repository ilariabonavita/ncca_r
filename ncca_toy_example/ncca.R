#' Nonparametric Canonical Correlation Analysis
#'
#' @description Compute the Nonparametric canonical correlations between two data matrices
#'
#' @importFrom FNN get.knnx()
#' @importFrom Matrix spMatrix(), Diagonal()
#' @importFrom irlba irlba()
#'
#' @param X The training data set containing n samples and p features
#' @param Y The paired training data set containing n samples and q features
#' @param XV The validation data set containing nv samples and p features
#' @param YV The paired validation data set containing nv samples and q features
#' @param XT The testing data set containing nt samples and p features
#' @param YT The paired testing data set containing nt samples and q features
#' @param d The dimension of the ouput transformed features
#' @param hx The bandwidth parameters for the KDEs of X
#' @param hy The bandwidth parameters for the KDEs of Y
#' @param nx The number of nearest neighbors for the KDEs of X
#' @param ny The number of nearest neighbors for the KDEs of Y
#' @param PreComputedNNs The name of the .Rda file containing two n*nx matrices for the nearest neighbor indice and the nearest neighbor Euclidean distances for X, and two n*ny matrices for the nearest neighbor indice and the nearest neighbor Euclidean distances for Y.  It is optional, only for debugging purposes.
#' @param ...
#'
#' @return X_new The d-dimensional projections of the training data X
#' @return Y_new The d-dimensional projections of the training data Y
#' @return XV_new The d-dimensional projections of the validation data X
#' @return YV_new The d-dimensional projections of the validation data Y
#' @return XT_new The d-dimensional projections of the testing data X
#' @return YT_new The d-dimensional projections of the testing data Y
#' @return The canonical correlation between X_new and Y_new
#' @return The canonical correlation between XV_new and YV_new
#'
#' @export
#' @author Meiwen Jia, \email{meiwen_jia@psych.mpg.de}
#' @references Michaeli, T., Wang, W., & Livescu, K. (2015). Nonparametric canonical correlation analysis. arXiv preprint arXiv:1511.04839.
#' @examples
#'
#' N <- 10000 # Overal number of examples (train+test)
#' N_paired <- 5000 # Number of training examples
#' MaxAngle <- 4*pi
#' MinRadius <-0.3
#' MaxRadius <- 8
#' NumNNs_X <- 20
#' NumNNs_Y <- 20
#' sx <- 0.5
#' sy <- 0.5
#' set.seed(8409)
#' ## Generate data for views 1,2
#' t <- seq(0, MaxAngle, length.out = N)
#' r <- seq(MinRadius, MaxRadius, length.out = N) + 2*runif(N)
#' #### generate X, the noise can be added!
#' X <- cbind(r*cos(t+0*rnorm(N)*0.05), r*sin(t+0*rnorm(N)*0.05))
#' X <- X + 0*matrix(rnorm(N*2), ncol = 2)
#' #### generate Y, the noise can be added!
#' Y <- cbind(t+0*rnorm(N)*1, 2*rnorm(N))
#' Y <- Y + 0*cbind(rep(0, N), rnorm(N))
#' ## Training data
#' PairedIndices <- sample(1:N, N_paired)
#' ## Test (or validation) data
#' UnpairedIndices <- setdiff(1:N,PairedIndices)
#' ncca_res <- ncca(X[PairedIndices,],Y[PairedIndices,], X[UnpairedIndices,],Y[UnpairedIndices,],
#'                  d = 2, hx = 0.75, hy = 0.75, nx = NumNNs_X, ny=NumNNs_Y)
#' cat("The nonparametric canonical correlation between X and Y is ", ncca_res$cor_XY, "\n")
#'
ncca <- function(X, # the paried training examples (rows) of features (columns)
                 Y, # the paried training examples (rows) of features (columns)
                 XV = NULL, # the paried validation examples
                 YV = NULL, # the paried validation examples
                 XT = NULL, # the paried test examples
                 YT = NULL, # the paried test examples
                 d, # the dimension of the output transformed features
                 hx = 1, # bandwith parameters for the KDEs of X
                 hy = 1, # bandwith parameters for the KDEs of Y
                 nx = 20, # number of nearest neighbours of X
                 ny = 20, # number of nearest neighbours of Y, it causes problem when nx != ny
                 PreComputedNNs = NULL, # only for debugging purposes
                 verbose = getOption("verbose"),
                 ...){

  ## Set default values to the unspecified parameters


  N <- nrow(X)
  ## Normalize data, using only training portion
  NXV <- nrow(XV); NXT <- nrow(XT);
  X <- rbind(X, XV, XT)
  NX <- nrow(X)
  rm(XV); rm(XT);

  NYV <- nrow(YV); NYT <- nrow(YT);
  Y <- rbind(Y, YV, YT)
  NY <- nrow(Y)
  rm(YV); rm(YT);

  X <- t(t(X) - colMeans(X[1:N, , drop = F]))
  X <- X/sqrt(mean(rowSums(X[1:N, , drop = F]^2)))

  Y <- t(t(Y) - colMeans(Y[1:N, , drop = F]))
  Y <- Y/sqrt(mean(rowSums(Y[1:N, , drop = F]^2)))

  ## compute NNs
  if(is.null(PreComputedNNs)){
    X_NNs <- FNN::get.knnx(X[1:N, ], X, k = nx)
    idxs_X <- t(X_NNs$nn.index)
    dists_X <- t(X_NNs$nn.dist)

    Y_NNs <- FNN::get.knnx(Y[1:N, ], Y, k = ny)
    idxs_Y <- t(Y_NNs$nn.index)
    dists_Y <- t(Y_NNs$nn.dist)
    rm(list = c("X_NNs","Y_NNs", "X", "Y"))
  }else{
    rm(list = c("X", "Y"))
    load(PreComputedNNs)
  }
  gc()

  ## compute weight matrices for X and Y
  colInd <- rep(c(1:NX), each = nx)
  Dx <- Matrix::spMatrix(nrow = N, ncol = NX, i = as.vector(idxs_X), j = colInd, x = exp(-0.5*as.numeric(dists_X)/(hx^2)))
  Dx <- t(Dx)
  colInd <- rep(c(1:NY), each = nx)
  Dy <- Matrix::spMatrix(nrow = N, ncol = NY, i = as.vector(idxs_Y), j = colInd, x = exp(-0.5*as.numeric(dists_Y)/(hx^2)))

  ## normalize the weight matrices
  Dx <- Matrix::Diagonal(NX, 1/rowSums(Dx))%*%Dx
  Dy <- Dy%*%t(Matrix::Diagonal(NY, 1/colSums(Dy)))

  ## Doubly stochastic normalization
  ## need some interpretation
  S <- Dx %*% Dy
  if(verbose) cat("Normalizing S to be doubly stochastic ...\n")
  its <- 15 # get error in SVD when the number of iterations increasing
  if(verbose) pb <- txtProgressBar(min = 0, max = its, style = 3)
  for(i in 1:its){
    S <- (NY/N) * Matrix::Diagonal(NY, 1/rowSums(S)) %*% S
    S <- S %*% Matrix::Diagonal(NX, 1/colSums(S)) * (NX/N)
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) close(pb)

  ## compute projections using the eigendecomposition of the kernel
  if(verbose) cat('Performing exact SVD ...\n')
  # S_SVD <- svd(S, d+1, d+1)
  S_SVD <- irlba::irlba(S, d+1) # a fast svd on a sparse matrix

  X_proj <- S_SVD$u[,2:(d+1)] * sqrt(NX)
  Y_proj <- S_SVD$v[,2:(d+1)] * sqrt(NY)

  ## output
  X_new <- X_proj[1:N, ]
  Y_new <- Y_proj[1:N, ]
  cor_XY <- cancor(X_new, Y_new)$cor

  if(is.null(NXV) | is.null(NYV)){
    XV_new <- YV_new <- cor_XVYV <- NULL
  }else{
    XV_new <- X_proj[N+(1:NXV), ]
    YV_new <- Y_proj[N+(1:NYV), ]
    cor_XVYV <- cancor(XV_new, YV_new)$cor
  }



  if(is.null(NXT) | is.null(NYT)){
    XT_new <- YT_new <- cor_XTYT <- NULL
  }else{
    XT_new <- X_proj[(N+NXV+1):NX, ]
    YT_new <- Y_proj[(N+NYV+1):NY, ]
    cor_XTYT <- cancor(XT_new, YT_new)$cor
  }

  ncca_out <- list(X_new = X_new, Y_new = Y_new,
                   XV_new = XV_new, YV_new = YV_new,
                   XT_new = XT_new, YT_new = YT_new,
                   cor_XY = cor_XY, cor_XVYV = cor_XVYV, cor_XTYT = cor_XTYT)
  return(ncca_out)
}

