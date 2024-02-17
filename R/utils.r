# local functions


wcor <- function(X, Y=X, w = rep(1,nrow(X))){
  # weighted correlation between matrix X and Y
  w <- w/sum(w)
  Xstd <- standardize_w(X, w)
  Ystd <- standardize_w(Y, w)
  t(Xstd) %*% diag(w) %*% Ystd
}

mean_w <- function(X,w = rep(1/nrow(X),nrow(X))){t(w/sum(w))%*% X}
center_w <- function(X,w = rep(1/nrow(X),nrow(X))){ X - rep(1,length(w))%*%t(w)%*% X }
standardize_w <- function(X,w = rep(1/nrow(X),nrow(X)), wsvd = FALSE){
  # NB requires w to be have sum 1
  ones <- rep(1,length(w))
  Xc <- X - ones %*% t(w)%*% X
  Xstd <- Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w))
  if (wsvd) Xstd <- wSVD(Xstd, w)
  return(Xstd)
}
