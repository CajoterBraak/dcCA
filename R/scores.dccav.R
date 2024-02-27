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
#' @param which_cor character or list of trait and environmental variables names (in this order)
#' in the data frames for which inter-set correlations must calculated.
#' Default: a character ("in_model") for all traits and variables in the model,
#' including collinear variables and levels.
#' @param scaling numeric (1,2 or 3) or character \code{"sites", "species" or "symmetric"}. Default: "sym".
#' Either site- (1) or species- (2) related scores are scaled by eigenvalues,
#' and the other set of scores is left unscaled,
#' or with 3 both are scaled symmetrically by square root of eigenvalues. Negative values are treated as the
#' corresponding positive ones by \code{abs(scaling)}. See also
#' @param tidy Return scores that are compatible with \code{ggplot2}:
#'  all scores are in a single data.frame, score type is identified by factor variable \code{score},
#'  the names by variable \code{label}, and species weights (in dc_CA_vegan) are in variable \code{weight}.
#'  See \code{\link[vegan]{scores.cca}}.
#' @param ...  Other arguments passed to the function (currently ignored).
#' @details
#'  An example of which_cor is: \code{which_cor = list(traits= c("SLA"), env = c("acidity","humidity") )}
#'  The function is modeled after \code{\link[vegan]{scores.cca}}.
#' @example demo/dune_dcCA.R
#' @returns A data frame if \code{tidy = TRUE}, a matrix if a single item is asked for and a named list of matrices if more than one item
#' is asked for. The following names can be included: \code{c("sites",
#' "constraints_sites", "centroids", "regression", "correlation", "biplot",
#' "species", "constraints_species", "regression_traits", "correlation_traits",
#' "biplot_traits", "centroids_traits")}. Each matrix has an attritute \code{"meaning"} explaining its content.
#' @export
scores.dccav <- function(x, choices=c(1,2), display= c("all"), scaling = "sym", which_cor = "in model", tidy = FALSE,...){
 # internal function
  f_meaning <- function( type_of_scores, scaling, txt){

    if (type_of_scores %in%
        c("sites", "constraints") )
      {point_type = "site"; opt_scal <- "sites"}
    else if (type_of_scores %in%  c("species", "constraints_species") )
      {point_type = "species"; opt_scal <- "species"}
    else if (type_of_scores %in% "centroids")
      {point_type = "environmental category"; opt_scal <- "sites"}
    else if (type_of_scores %in% "centroids_traits")
    {point_type = "trait category"; opt_scal <- "species"}
    else if (type_of_scores %in% c("biplot", "biplot_traits"))
    {point_type = "arrows"; opt_scal <- scaling}

    txt1 <-   paste(txt," in scaling '", scaling, "' optimal for biplots and inter-", point_type, " distances.", sep = "")
    txt1b <-  paste(txt," in scaling '", scaling, "' optimal for biplots and, almost so, for inter-", point_type, " distances.", sep = "")
    txt1c <-  paste(txt," in scaling '", scaling, "' optimal for biplots, but unsuited for inter-", point_type, " distances.", sep = "")

    txt2<-    paste(txt," in scaling '", scaling, "', optimal for biplots displays, suboptimal for display of inter-", point_type, " distances.", sep = "")
    txt3 <-   paste(txt," in scaling '", scaling, "' optimal for biplot displays, suboptimal for distance interpretation.", sep = "")
    txt4 <-   paste(txt," in scaling '", scaling, "'.", sep = "")

    if (type_of_scores %in% c("biplot", "biplot_traits")) txt1 <- txt2<- txt3 <-txt4


    thres1 <- 1.5 # to be decided upon... What does Canoco 5 suggest?
    thres2 <- 4
    ratio_eig <- x$eigenvalues[choices[1]]/x$eigenvalues[choices[2]]
    if (scaling == "symmetric"){
      if (sqrt(ratio_eig) < thres1  ) txt3 <- txt1b else if (sqrt(ratio_eig) > thres2) txt3 <- txt1c
    } else
      if (ratio_eig < thres1 ) txt2 <- txt1b else if (ratio_eig > thres2) txt2 <- txt1c


    if (opt_scal == scaling) txt_out <- txt1
    else if (scaling == "symmetric") txt_out <- txt3
    else txt_out <- txt2

    return(txt_out)

  }
 #
 if (!class(x)[1]=="dccav") stop("The first argument must be the result of the function dc_CA_vegan.")

 tabula <- c( "sites", "constraints", "regression", "biplot", "correlation",
         "centroids","species", "constraints_species", "regression_traits","biplot_traits" ,
             "correlation_traits","centroids_traits" )
 names(tabula) <- c( "wa", "lc", "reg","bp", "cor", "cn","sp","lc_traits", "reg_traits","bp_traits", "cor_traits","cn_traits")
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
  # make sure axes chosen by choices are not larger than the rank
  choices <- choices[choices <= length(x$eigenvalues)]

  if (is.character (scaling) ){
    scaling <- match.arg(scaling,  c( "sites", "species","symmetric"))
    if (scaling == "sites") num_scaling <- 1 else if (scaling == "species") num_scaling <- 2 else if (scaling =="symmetric") num_scaling = 3 else stop("scaling type not recognized")
  } else if (is.numeric(scaling)) {
      num_scaling <- abs(scaling)
      scaling <- c("sites","species","symmetric")[num_scaling]
  } else stop("scaling type not recognized")

  slam <- sqrt(x$eigenvalues[choices])
  scal <- list(rep(1, length(slam) ), slam, sqrt(slam))[[abs(num_scaling)]]
  diag_scal_sites   <- diag(1/scal)
  diag_scal_species <- diag(scal)



    sol <- list()


# scaling for site related scores (incl env) ------------------------------



    if ( "sites" %in%take){
      sol$sites <-   site_axes$site_scores[[1]][,choices, drop = FALSE] %*% diag_scal_sites
      attr(sol$sites, which = "meaning") <-f_meaning("sites", scaling,
      "CMWs of the trait axes (constraints species)")
    }

    if ( "constraints" %in%take){
      sol$constraints_sites <- site_axes$site_scores[[2]][,choices, drop = FALSE] %*% diag_scal_sites
      attr(sol$constraints_sites, which = "meaning") <- f_meaning("constraints", scaling,
      paste("linear combination of the environmental predictors",
        "and the covariates (making the ordination axes orthogonal to the covariates)'", collapse = "")
      )
    }

    if ( "centroids" %in%take){

      if (which_cor == "in model") {
        in_model <- get_focal_and_conditioning_factors(x$RDAonEnv, factors_only = FALSE)$`focal factor`
      } else in_model = which_cor
      dat = x$data$dataEnv[, in_model, drop= FALSE]
      cn  <- centroids.cca(x$site_axes$site_scores$site_scores_unconstrained,
                           dat, wt=x$weights$rows)[,choices, drop = FALSE]

      if(!is.null(cn)){
        cn <- cn %*% diag_scal_sites
        attr(cn, which = "meaning") <- f_meaning("centroids", scaling,
                                                 "environmental category means of the site scores")
      }
      sol$centroids <-  cn
    }

    if ("regression"%in% take) {

      regr <- c_env_normed[,choices +3] %*% diag_scal_sites
      if (tidy)sol$regression <- regr else
      sol$regression <-cbind(c_env_normed[,1:3], regr)

      attr(sol$regression, which = "meaning")<-
        paste("mean, sd, VIF, standardized regression coefficients and their optimistic t-ratio in scaling '",scaling,"'.",sep="")
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
      sol$biplot <- e_rcor%*% diag(slam/R) %*% diag_scal_sites
      colnames(sol$biplot)<- paste("dcCA", choices, sep = "")
      attr(sol$biplot, which = "meaning") <- f_meaning("biplot", scaling,
      "biplot scores of environmental variables for display with biplot-traits for fourth-corner correlations")
    }






# Species stats -----------------------------------------------------------
   if ( "species" %in%take) {
     sol$species <- species_axes$species_scores[[1]][,choices, drop = FALSE] %*% diag_scal_species
     attr(sol$species, which = "meaning")<-f_meaning("species", scaling,
     "SNC on the environmental axes (constraints sites)")

   }
   if ( "constraints_species" %in%take){
     sol$constraints_species <- species_axes$species_scores[[2]][,choices, drop = FALSE] %*% diag_scal_species
     attr(sol$constraints_species, which = "meaning")<- f_meaning("constraints_species", scaling,
       paste("linear combination of the traits",
      "and the trait covariates (making the ordination axes orthogonal to the covariates)'", collapse = "")
     )

   }
    if ("regression_traits"%in% take){
      regr <- species_axes$c_traits_normed[,choices + 3] %*% diag_scal_species
      if (tidy) sol$regression_traits <-  regr else
      sol$regression_traits <- cbind(species_axes$c_traits_normed[, 1:3] , regr)
      attr(sol$regression_traits, which = "meaning")<-
        paste("mean, sd, VIF, standardized regression coefficients and their optimistic t-ratio in scaling '",scaling,"'.",sep="")

    }
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
      sol$biplot_traits <- t_rcor %*% diag(diag(diag_scal_species)/R)
      colnames(sol$biplot_traits)<- paste("dcCA", choices, sep = "")
      attr(sol$biplot_traits, which = "meaning") <-f_meaning("biplot", scaling,
        "biplot scores of traits for display with biplot scores for fourth-corner correlation")
    }

    if ( "centroids_traits" %in%take){

      if (which_cor == "in model") {
        # in_model <- colnames(x$data$dataTraits)%in% rownames(attr(stats::terms(x$CCAonTraits), which = "factors"))
        in_model <- get_focal_and_conditioning_factors(x$CCAonTraits, factors_only = FALSE)$`focal factor`
      } else in_model = which_cor
      dat = x$data$dataTraits[, in_model, drop= FALSE]
      cn  <- centroids.cca(x$species_axes$species_scores$species_scores_unconstrained,
                      dat, wt=x$weights$columns)[,choices, drop = FALSE]

      # sol$centroids_traits_lc <-  centroids.cca(x$species_axes$species_scores$lc_traits_scores,
      #                                       dat, wt=x$weights$columns)[,choices, drop = FALSE] %*% diag_scal_species
      if(!is.null(cn)){
        cn <- cn %*% diag_scal_species
        attr(cn, which = "meaning") <-  f_meaning("centroids_traits", scaling,
                                                  "trait category means of the species scores")
      }

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
     attr(sol, which = "scaling") <- scaling
   }


    ## return NULL instead of list(), and matrix instead of a list of
    ## one matrix
    switch(min(2, length(sol)), sol[[1]], sol)
}
