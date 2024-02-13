
data("dune_trait_env")

rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
out <- dc_CA_vegan(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~. ,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv =dune_trait_env$envir,
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
                   dataTraits =dune_trait_env$traits[,-c(1,2)],
                   verbose = TRUE)
mod_scores <- scores(out, display = "all") # not yet exported
names(mod_scores)
#[1] "sites"      "lc"         "species"    "lc_traits"  "cor"
#[6] "reg_traits" "cor_traits"

# community-level permutation test
anova(out$RDAonEnv) # all option of anova.cca are available!
# a species-level permuation test requires an dedicated new function
# this is not yet available in this version

# for illustration: a dc-CA model with a trait covariate
out2 <- dc_CA_vegan(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~ SLA+Height+ LDMC+ Lifespan +Condition(Seedmass) ,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv =dune_trait_env$envir,
                   dataTraits =dune_trait_env$traits,
                   verbose = TRUE)

# for illustration: a dc-CA model with both environmental and trait covariates
out3 <- dc_CA_vegan(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                    formulaTraits = ~ SLA+Height+LDMC+Lifespan +Condition(Seedmass) ,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv =dune_trait_env$envir,
                    dataTraits =dune_trait_env$traits,
                    verbose = TRUE)

# for illustration: same model but using out2 for speed, as the trait model and data did not change
out3B <- dc_CA_vegan(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                    dataEnv =dune_trait_env$envir,
                    dc_CA_vegan_object = out2,
                    verbose = TRUE)
all.equal(out3,out3B) # TRUE

# All statistics and scores have been checked against the results with "focus on Case distances" (=Sites)
# in Canoco 5.15 (ter Braak & Smilauer, 1918).

