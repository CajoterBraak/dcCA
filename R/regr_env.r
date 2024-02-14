#' Standardized regression coefficient with optimistic t-ratios and
#' mean and standard deviation of the environmental variables
#'
#' @noRd
regr_env <-function(out){

if (!class(out)[1]=="dccav") stop("The first argument must be the result of the function dc_CA_vegan.")

step2 <- out$RDAonEnv
c_env_normed <- as.matrix(vegan:::scores.rda(step2, display = "reg", scaling = "species", choices = 1:step2$CCA$rank))
c_env_normed <-   c_env_normed %*% diag(sqrt(step2$CCA$eig)) # in sites scaling
colnames(c_env_normed) <- paste("Regr", seq_len(ncol(c_env_normed)), sep = "")
avgX <- attr(step2$CCA$QR$qr, which= "scaled:center")
sdsX = sqrt(colMeans(qr.X(step2$CCA$QR)^2))

## t-values of regression coefficients based on type = "canoco" residuals
tval <- coef(step2)/sqrt(diag(vegan:::vcov.cca(step2, type = "canoco")))

VIF <- fVIF(step2$CCA$QR) #diag(XtX_inv)*sdsX^2

colnames(tval) <- paste("tval", seq_len(ncol(tval)), sep = "")
c_env_normed <- cbind(Avg = avgX, SDS = sdsX, VIF = VIF, c_env_normed[names(sdsX),],tval[names(sdsX),])
attr(c_env_normed, "meaning") <- "mean, sd, VIF, standardized regression coefficients and their optimistic t-ratio"

return(c_env_normed)
}

fVIF <- function(object) {
  # object is a qr
  if(is.qr(object)){
    R <- qr.R(object)
    R_inv <- solve(R)
    XtX_inv <- R_inv %*% t(R_inv)
    VIF <- diag(XtX_inv)* colSums(qr.X(object)^2)
    } else VIF <- NA
  return(VIF)
}
