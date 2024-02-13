#' @title Extract results of a vegan-based double constrained correspondence Analysis (dc-CA)
#'
#' @description
#' This function works very much like the \code{vegan} \code{\link[vegan]{scores}} function,
#' in particuluar \code{\link[vegan]{scores.cca}}, with the additional results such
#' as regression coefficients and linear combinations of traits ('regr_traits','lc_traits').
#' In the current version, there is a single scaling (scaling = "sites").
#' All scores from CCA obey the so called transition formulas and so do the scores of dc-CA.
#' The difference is that the linear combinations of traits (the \emph{constrained} species scores)
#' replace the usual (\emph{unconstrained}) species scores.
#'
#' @param x object
#' @param choices integer vector of which axes to obtain. Default: all dc-CA axes.
#' @param display a character vector, one or more of
#' \code{c("all","species","sites","sp", "wa", "lc", "cor",
#' "reg", "cn","lc_traits", "reg_traits", "cor_traits")}.
#' The first nine are as in \code{\link[vegan]{scores.cca}} and remaining ones are similar scores for traits.
#' @param which_cor character or list of trait and environmental variables names in the data frames
#' for which inter-set correlations must calculated.
#' Default: a character ("in_model") for all traits and variables in the model,
#' including collinear variables and levels.
#'
#' @example demo/dune_dcCA.R
#' @export
scores.dccav <- function(x, choices, display= c("all"), which_cor = "in model", tidy = FALSE,...){
 scaling = "sites" # currently only a single scaling is available
 if (!class(x)[1]=="dccav") stop("The first argument must be the result of the function dc_CA_vegan.")

 tabula <- c("species", "sites", "constraints", "correlation",
             "regression", "centroids", "constraints_traits", "regression_traits", "correlation_traits" )
 names(tabula) <- c("sp", "wa", "lc", "cor", "reg", "cn","lc_traits", "reg_traits", "cor_traits")
  display <- match.arg(display, c("sites", "species", "wa",
                                 "lc", "cor", "reg", "cn","lc_traits", "reg_traits", "cor_traits", "all"),
                      several.ok = TRUE)
 ## set "all" for tidy scores
 if (tidy)
   display <- "all"
 if("sites" %in% display)
   display[display == "sites"] <- "wa"
 if("species" %in% display)
   display[display == "species"] <- "sp"
  if("correlation" %in% display)
    display[display == "correlation"] <- "co"
 if("all" %in% display)
   display <- names(tabula)
 take <- tabula[display]

  if ((!"species_axes"%in%names(x)) && any(c("species",
           "constraints_traits", "regression_traits", "correlation_traits")%in% take)){
    c_env_normed <- regr_env(out)
    species_axes <- f_trait_axes(out)
  } else if ("species_axes"%in%names(x)){c_env_normed <- x$c_env_normed; species_axes<- x$species_axes}

  if (scaling == "sites")  myconst <- sqrt(vegan:::nobs.cca(x$RDAonEnv)*x$RDAonEnv$tot.chi) else
    if (scaling == "species") myconst <- sqrt(vegan:::nobs.cca(x$RDAonEnv))

   sol <- list()

    if ("sites" %in% take)
    sol$sites  <- vegan:::scores.rda(x$RDAonEnv, display = c("sites"), scaling = scaling,
                                       choices = seq_len(x$RDAonEnv$CCA$rank), const = myconst)
    if ( "constraints" %in%take)
    sol$lc  <- vegan:::scores.rda(x$RDAonEnv, display = c("lc"), scaling = scaling,
                                    choices = seq_len(x$RDAonEnv$CCA$rank), const = myconst)



    if ( "species" %in%take) sol$species <- species_axes$species[[1]]
    if ( "constraints_traits" %in%take)sol$lc_traits <- species_axes$species[[2]]

    if ("regression"%in% take) sol$reg <- c_env_normed
    if ("correlation"%in% take) {
      sites  <- vegan:::scores.rda(x$RDAonEnv, display = c("sites"), scaling = scaling,
                                   choices = seq_len(x$RDAonEnv$CCA$rank))
      # correlations of the dataEnv wrt the first axis (site scores)
      if (!is.list(which_cor )) in_model <- colnames(x$data$dataEnv)%in% colnames(attr(terms(x$RDAonEnv), which = "factors")) else
        in_model = which_cor[[2]]
      env0 <-  model.matrix(~.-1, constrasts = FALSE, data = x$data$dataEnv[, in_model, drop = FALSE])
      Cormat <- cov2cor(cov(cbind( env0, sites)))
      Cor_Env_CWM <- Cormat[seq_len(ncol(env0)),ncol(env0) + seq_len(x$RDAonEnv$CCA$rank) , drop = FALSE]
      colnames(Cor_Env_CWM) <- paste("CWM-ax", seq_len(ncol(Cor_Env_CWM)), sep= "")
      attr(Cor_Env_CWM, which = "meaning")<-
"inter set correlation, correlation between environmental variables and the sites scores (CWMs)"
      sol$cor <- Cor_Env_CWM
    }

    if ("regression_traits"%in% take)sol$reg_traits <- species_axes$c_traits_normed
    if ("correlation_traits"%in% take) {
      if (!is.list(which_cor)){ sol$cor_traits <- species_axes$correlation} else{
        whichc = which_cor[[1]]
        cor_traits_SNC <- f_trait_axes(out, which_cor = whichc)
        sol$cor_traits <- cor_traits_SNC$correlation
      }
    }



    ## return NULL instead of list(), and matrix instead of a list of
    ## one matrix
    switch(min(2, length(sol)), sol[[1]], sol)
}
