#' Obtain Column (Species) Related Scores
#' @param out object from \code{\link{dc_CA_vegan}}
#' @param which_cor character or names of traits
#' for which inter-set correlations must calculated.
#' Default: "in_model" for all traits and variables in the model specified by \code{formulaTraits}.
#' @noRd
# @export
f_trait_axes <- function(out, which_cor = "in model"){
  # which_cor character or names of traits
  # for which inter-set correlations must calculated.
  #  Default "in_model" for all traits and variables in the model
  # SNC lc_traits and trait regr, tval, cor
  N <- stats::nobs(out$RDAonEnv)
  lc_scores  <- vegan::scores(out$RDAonEnv, display = c("lc"), scaling = "species",
                                   choices = seq_len(Rank_mod(out$RDAonEnv)), const = sqrt(N))

  SNC <-  (t(as.matrix(out$data$Y)) %*% lc_scores) / (out$CCAonTraits$rowsum * N)

  if (!is.null(out$CCAonTraits$pCCA)){  # orthogalize with respect to any covariate
    SNC <- SNC - calculate_b_se_tval(get_QR(out$CCAonTraits, model= "pCCA"), y=SNC,
                                     w = stats::weights(out$CCAonTraits, "sites"),  scale2 = 0, name = "SNC", fitted_only = TRUE)
  }
  #print(names(out))
  res <- calculate_b_se_tval(get_QR(out$CCAonTraits), y=SNC,
                             w = stats::weights(out$CCAonTraits, "sites"),  scale2 = 1, name = "SNC")

  #print(names(res))

  c_traits_normed <- res$coef_normed
  attr(c_traits_normed, which = "warning") <-"The t-values are optimistic, i.e. an underestimate of their true absolute value"


  # correlations of the dataTraits with the SNC wrt the axes
  if (which_cor == "in model") {
    in_model <- colnames(out$data$dataTraits)%in% rownames(attr(stats::terms(out$CCAonTraits), which = "factors"))

    } else in_model = which_cor
  traits0 <-  stats::model.matrix(~.-1, constrasts = FALSE, data = out$data$dataTraits[, in_model, drop= FALSE])
  #traits0 <-  model.matrix(~. -1, constrasts = FALSE, data = out$data$dataTraits)
  Cormat <- stats::cov2cor(ade4::covwt(cbind( traits0, SNC), w= out$CCAonTraits$rowsum))
  Cor_Trait_SNC <- Cormat[seq_len(ncol(traits0)),ncol(traits0) + seq_len(Rank_mod(out$RDAonEnv)) , drop = FALSE]
  colnames(Cor_Trait_SNC) <- paste("SNC-ax", seq_len(ncol(Cor_Trait_SNC)), sep= "")
  attr(Cor_Trait_SNC, which = "meaning")<- "inter set correlation, correlation between traits and SNC of axes"

  out2 <- list(species_scores = list(species_scores_unconstrained = res$y,
                                     lc_traits_scores = res$fitted), correlation = Cor_Trait_SNC, c_traits_normed= c_traits_normed, b_se= res$b_se, R2_traits = res$R2)
  return(out2)
}
