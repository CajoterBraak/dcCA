
data("dune_trait_env")
# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
# must delete "Sites" from response matrix or data frame
Y <- dune_trait_env$comm[,-1] # must delete "Sites"

out <- dc_CA_vegan(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan ,
                   response = Y,
                   dataEnv =dune_trait_env$envir,
                   dataTraits =dune_trait_env$traits,
                   verbose = FALSE)
out <- print(out) # more efficient for  scores() than just 'print(out)' when verbose = FALSE

mod_scores <- scores(out, display = "all")
names(mod_scores)
#str(mod_scores)
axis <- 1
stats <- c(community= "cor", species = "cor") # "weights
if (length(stats)==1)stats <- c(stats,stats)


# community level plot ----------------------------------------------------

if (stats[1] == "weights"){
  trait_scores <- mod_scores$regression_traits
  ylab_traits <-  "Weight in composite trait"
  offset  <- 3
} else {
  trait_scores <- mod_scores$correlation_traits
  ylab_traits <- "cor(trait, SNC)"
  offset <- 0
}

#smoothing_method <- c(sites = "lm",species = "gam")
#if (length(smoothing_method)==1)smoothing_method <- c(smoothing_method,smoothing_method)

dat <- out$data$dataEnv
dat$composite_env <- mod_scores$constraints_sites[,axis]
dat$CWM_composite_trait <- mod_scores$sites[,axis]

#site_groups<- PRC::get_focal_and_conditioning_factors(out$RDAonEnv)$condition # site groups
site_groups<- NULL  # no site groups
if (!is.null(site_groups)) dat$site_groups <- dat[[site_groups]]

library(ggplot2)
if (is.null(site_groups))   p <- ggplot(data=  dat, aes(x = composite_env, y= CWM_composite_trait)) else
  p <- ggplot(data=  dat, aes(x = composite_env, y= CWM_composite_trait, group = .data[[site_groups]], color= .data[[site_groups]]))
p_env <- p + geom_point() + geom_smooth(method = lm) + xlab("environmental gradient")
# or to show that PB has no relevance
#pp <- p + geom_point() + geom_smooth()
#pp


plot_traits <- plot_species_scores_bk(
  species_scores= trait_scores,
  ylab = ylab_traits,
  threshold = 0,
  y_lab_interval = 0.2,
  speciesname = NULL,
  scoresname = colnames(trait_scores)[offset+ axis],
  selectname = "Fratio1",
  verbose = TRUE
)

# species-level plot ------------------------------------------------------

if (stats[2]=="weights"){
  env_scores <- mod_scores$regression
  ylab_env <-  "Weight in gradient"
  offset <- 3
} else {
  env_scores <- mod_scores$correlation
  ylab_env <- "cor(env, CWM)"
  offset <- 0
}


dat <- out$data$dataTraits
dat$composite_trait <- mod_scores$constraints_species[,axis]
dat$SNC_composite_env <- mod_scores$species[,axis]
dat$`rel.abundance` <- out$CCAonTraits$rowsum

#fc<- PRC::get_focal_and_conditioning_factors(out$RDAonEnv)
species_groups <- NULL # species groups
if(!is.null(species_groups)) dat$species_groups <- species_groups

library(ggplot2)
if (is.null( species_groups))   p <- ggplot(data=  dat, aes(x = composite_trait, y= SNC_composite_env )) else
  p <- ggplot(data=  dat, aes(x = composite_trait, y= SNC_composite_env, group = species_groups, color= species_groups, size=tot_abun))
#p_traits <- p + geom_point() + stat_smooth(method = lm) + xlab("trait composite") + ylab("SNC along gradient")

p_traits <- p + geom_point(aes( size  = `rel.abundance`)) +
  stat_smooth(aes(weight = `rel.abundance`),method = "lm") +
  xlab("trait composite") + ylab("SNC along gradient")+
  theme(legend.position  = c(0.1,0.75))

#p_traits
#suppressWarnings(print(p_traits))

plot_env <- plot_species_scores_bk(
  species_scores= env_scores,
  ylab = ylab_env,
  threshold = 0,
  y_lab_interval = 0.2,
  speciesname = NULL,
  scoresname = colnames(env_scores)[offset+ axis],
  selectname = "Fratio1",
  verbose = TRUE
)

# modifying the plot
#gg <- plotPRC(mod_prc, plot = "ditch",width = c(4,1), verbose = FALSE)
#p1 <- gg$separateplots$treatments + ggplot2::ggtitle(paste("new title:", latex2exp::TeX("$c_{dt}$")))   # PRC plot of samples (c_dt)
#p2 <- gg$separateplots$species    + ggplot2::ylab("new title: loadings")# loadings of species  (b_k)
# Assign these plots to symbols and use grid.arrange to produce the plot  you like, for example:

# plot arrange ------------------------------------------------------------


if (stats[1]=="weights"){
gridExtra::grid.arrange(p_env+ ylab("CWM of composite trait"),
                        plot_traits, ncol =2, widths = c(4,1),
                        top = "Community level of double constrained correspondence analysis",
                        left ="", right =  "")
suppressWarnings( print(gridExtra::grid.arrange(p_traits,
                        plot_env, ncol =2, widths = c(4,1),
                        top = "Species level of double constrained correspondence analysis",
                        left ="", right =  "")))

} else {
  gridExtra::grid.arrange(p_env+ ylab("CWM of composite trait"),
                          plot_env, ncol =2, widths = c(4,1),
                          top = "Community level of double constrained correspondence analysis",
                          left ="", right =  "")
  suppressWarnings( print(gridExtra::grid.arrange(p_traits,
                                                  plot_traits, ncol =2, widths = c(4,1),
                                                  top = "Species level of double constrained correspondence analysis",
                                                  left ="", right =  "")))

}






