#' @title Extract results of a vegan-based double constrained correspondence Analysis (dc-CA)
#'
#' @description
#' This function works very much like the \code{vegan} \code{\link[vegan]{scores}} function,
#' in particular \code{\link[vegan]{scores.cca}}, with the additional results such
#' as regression coefficients and linear combinations of traits \code{('regr_traits','lc_traits')}
#' In the current version, there is a single scaling (\code{scaling = "sites"}).
#' All scores from CA obey the so called transition formulas and so do the scores of CCA and dc-CA.
#' The differences are, for CCA, that the linear combinations of environmental variables
#' (the \emph{constrained} site scores)
#' replace the usual (\emph{unconstrained}) site scores, and for dc-CA,
#' that the linear combinations of traits (the \emph{constrained} species scores)
#' also replace the usual (\emph{unconstrained}) species scores.
#'
#' @param x object
#' @param choices integer vector of which axes to obtain. Default: all dc-CA axes.
#' @param display a character vector, one or more of
#' \code{c("all","species","sites","sp", "wa", "lc","bp", "cor", "reg", "cn",
#' "lc_traits", "reg_traits", "cor_traits","bp_traits","cn_traits")}.
#' The first ten are as in \code{\link[vegan]{scores.cca}} (except \code{"cor"})
#' and remaining ones are similar scores for traits.
#' @param which_cor character or list of trait and environmental variables names in the data frames
#' for which inter-set correlations must calculated.
#' Default: a character ("in_model") for all traits and variables in the model,
#' including collinear variables and levels.
#' @param tidy Return scores that are compatible with \code{ggplot2}:
#'  all scores are in a single data.frame, score type is identified by factor variable \code{score},
#'  the names by variable \code{label}, and species weights (in dc_CA_vegan) are in variable \code{weight}.
#'  See \code{\link[vegan]{scores.cca}}.
#' @param ...  Other arguments passed to the function (currently ignored).
# @details
#' @example demo/dune_dcCA.R
#' @export
scores.dccav <- function(x, choices=c(1,2), display= c("all"), which_cor = "in model", tidy = FALSE,...){
 #
 scaling = "sites" # currently only a single scaling is available
 #
 if (!class(x)[1]=="dccav") stop("The first argument must be the result of the function dc_CA_vegan.")

 tabula <- c("species", "sites", "constraints", "regression", "biplot", "correlation",
             "centroids", "constraints_species", "regression_traits","biplot_traits" ,
             "correlation_traits","centroids_traits" )
 names(tabula) <- c("sp", "wa", "lc", "reg","bp", "cor", "cn","lc_traits", "reg_traits","bp_traits", "cor_traits","cn_traits")
 #print("here is scores.dccav")
 display <- match.arg(display,
                      c("sp", "wa", "lc","bp", "cor", "reg", "cn","lc_traits", "reg_traits","bp_traits", "cor_traits","cn_traits","sites", "species", "all"),
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

  if (!"species_axes"%in%names(x)){
    site_axes <- f_env_axes(x)
    c_env_normed <- site_axes$c_env_normed
    species_axes <- f_trait_axes(x)
  } else if ("species_axes"%in%names(x)){
    c_env_normed <- x$c_env_normed; species_axes<- x$species_axes; site_axes<- x$site_axes
    }

  if (scaling == "sites")  myconst <- sqrt(x$Nobs*x$RDAonEnv$tot.chi) else
    if (scaling == "species") myconst <- sqrt(x$Nobs)
  # make sure axes chosen by choices are not larger than the rank
  choices <- choices[choices <= Rank_mod(x)]
  if (tidy) regchoices <-  choices+3 else regchoices <- c(1:3, choices+3) # coefs only (tidy) or with mean,sd,vif

    sol <- list()

    if ("sites" %in% take){
    # sol$sites  <- vegan::scores(x$RDAonEnv, display = c("sites"), scaling = scaling,
    #                                    choices = choices, const = myconst)
    sol$sites <- site_axes$site_scores[[1]][,choices, drop = FALSE]

    attr(sol$sites, which = "meaning") <- "CMWs of the trait axes (constraints species) in 'Sites' scaling."
    }
    if ( "constraints" %in%take){
    #sol$constraints_sites  <- vegan::scores(x$RDAonEnv, display = c("lc"), scaling = scaling,
    #                                choices = choices, const = myconst)
    sol$constraints_sites <- site_axes$site_scores[[2]][,choices, drop = FALSE]

    attr(sol$constraints_sites, which = "meaning") <- c("linear combination of the environmental predictors",
      "(and the covariates, so as to make the ordination axes orthogonal to the covariates)")
    }


   if ( "centroids" %in%take){
      cn <-  vegan::scores(x$RDAonEnv, display = c("cn"), scaling = scaling,
                                       choices = choices, const = myconst)
     if(!is.null(cn)){
     attr(cn, which = "meaning") <- "environmental category means of the ordination axes  (constraints sites)"

     }
    sol$centroids <- cn
   }


    if ("regression"%in% take) {
      sol$regression <- c_env_normed[,regchoices]
      attr(sol$regression, which = "meaning")<-
        "mean, sd, VIF, standardized regression coefficients and their optimistic t-ratio"
    }
    if ("correlation"%in% take) {
        if (!is.list(which_cor)){
          sol$correlation <- site_axes$correlation[,choices, drop = FALSE]
        } else{
          whichc = which_cor[[2]]
          cor_Env_CWM <- f_env_axes(x, which_cor = whichc)
          sol$correlation <- cor_Env_CWM$correlation[,choices, drop = FALSE]
        }
      attr(sol$correlation,  which = "meaning")<-
        "inter set correlation, correlation between environmental variables and the sites scores (CWMs)"
    }

    if ( "biplot" %in%take){
      e_rcor <- x$site_axes$correlation[,choices, drop = FALSE]
      R <- sqrt(x$site_axes$R2_env[choices])
      sing <- sqrt(x$eigenvalues[choices])
      sol$biplot <- e_rcor%*% diag(sing/R)
      colnames(sol$biplot)<- paste("dcCA", choices, sep = "")
      attr(sol$biplot, which = "meaning") <-
  "biplot scores of environmental variables for display with biplot-traits for fourth-corner correlations"
    }






# Species stats -----------------------------------------------------------
   if ( "species" %in%take) {
     sol$species <- species_axes$species_scores[[1]][,choices, drop = FALSE]
     attr(sol$species, which = "meaning")<- "SNC on the ordination axes (constraints sites), scaled to unit weighted sum of squares"
   }
   if ( "constraints_species" %in%take){
     sol$constraints_species <- species_axes$species_scores[[2]][,choices, drop = FALSE]
     attr(sol$constraints_species, which = "meaning")<- c("linear combination of the traits",
        "(and the trait covariates, so as to make the ordination axes orthogonal to the trait covariates)")
   }
    if ("regression_traits"%in% take)sol$regression_traits <- species_axes$c_traits_normed[,regchoices]
    if ("correlation_traits"%in% take) {
      if (!is.list(which_cor)){
        sol$correlation_traits <- species_axes$correlation[,choices, drop = FALSE]
      } else{
        whichc = which_cor[[1]]
        Cor_Trait_SNC <- f_trait_axes(x, which_cor = whichc)
        sol$correlation_traits <- Cor_Trait_SNC$correlation[,choices, drop = FALSE]
      }
      attr(sol$correlation_traits, which = "meaning")<-
        "inter set correlation, correlation between traits and the species scores (SNCs)"
    }

    if ( "biplot_traits" %in%take){
      t_rcor <- x$species_axes$correlation[,choices, drop = FALSE]
      R <- sqrt(x$species_axes$R2_traits[choices])
      sol$biplot_traits <- t_rcor %*% diag(1/R)

      colnames(sol$biplot_traits)<- paste("dcCA", choices, sep = "")
      attr(sol$biplot_traits, which = "meaning") <-
        "biplot scores of traits for display with biplot scores for fourth-corner correlations"
    }

    if ( "centroids_traits" %in%take){

      if (which_cor == "in model") {
        # in_model <- colnames(x$data$dataTraits)%in% rownames(attr(stats::terms(x$CCAonTraits), which = "factors"))
        in_model <- get_focal_and_conditioning_factors(x$CCAonTraits)$`focal factor`
      } else in_model = which_cor
      dat = x$data$dataTraits[, in_model, drop= FALSE]
      cn  <- centroids.cca(x$species_axes$species_scores$species_scores_unconstrained,
                      dat, wt=x$weights$columns)[,choices, drop = FALSE]
      sol$centroids_traits_lc <-  centroids.cca(x$species_axes$species_scores$lc_traits_scores,
                                            dat, wt=x$weights$columns)[,choices, drop = FALSE]
      if(!is.null(cn))
        attr(cn, which = "meaning") <- "trait category means of the ordination axes  (constraints sites)"
      sol$centroids_traits <-  cn
    }

   for (nam in names(sol)){
     if(!is.null(sol[[nam]])){
     if (!nam %in% c("regression", "regression_traits", "correlation", "correlation_traits"))
     colnames(sol[[nam]]) <-  paste("dcCA", choices, sep = "") else if (nam %in% c("regression", "regression_traits"))
       colnames(sol[[nam]])[-c(1,2,3)] <-  paste("dcCA", choices, sep = "")
     }
   }

# taken from vegan::scores.cca with thanks --------------------------------
   ## Take care that scores have names
   if (length(sol)) {
     for (i in seq_along(sol)) {
       if (is.matrix(sol[[i]]))
         rownames(sol[[i]]) <-
           rownames(sol[[i]], do.NULL = FALSE,
                    prefix = substr(names(sol)[i], 1, 3))
     }
   }
   ## tidy scores
   if (tidy) {
     if (length(sol) == 0) # no requested scores existed
       return(NULL)
     ## re-group biplot arrays duplicating factor centroids
     if (!is.null(sol$biplot) && !is.null(sol$centroids)) {
       dup <- rownames(sol$biplot) %in% rownames(sol$centroids)
       if (any(dup)) {
         sol$factorbiplot <- sol$biplot[dup,, drop=FALSE]
         sol$biplot <- sol$biplot[!dup,, drop=FALSE]
       }
     }
     group <- sapply(sol, nrow)
     group <- rep(names(group), group)
     sol <- do.call(rbind, sol)
     label <- rownames(sol)
     cw <- x$weights$columns
     rw <- x$weights$rows
     w <- rep(NA, nrow(sol))
     if (any(weighted <- group == "sites"))
       w[weighted] <- rw
     if (any(weighted <- group == "constraints"))
       w[weighted] <- rw
     if (any(weighted <- group == "species"))
       w[weighted] <- cw
     if (any(weighted <- group == "constraints_species"))
       w[weighted] <- cw
     sol <- as.data.frame(sol)
     sol$score <- as.factor(group)
     sol$label <- label
     sol$weight <- w
     names(sol)[seq_along(choices)] <- paste("dcCA", choices, sep = "")
   }


    ## return NULL instead of list(), and matrix instead of a list of
    ## one matrix
    switch(min(2, length(sol)), sol[[1]], sol)
}
