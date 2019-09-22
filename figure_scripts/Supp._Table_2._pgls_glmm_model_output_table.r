#Table of pgls and lm parameter estimates and 95% credible intervals.
rm(list=ls())
source('paths.r')
library(expss)

#load data.---
pg <- readRDS(pgls.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.path)
pg <- pg[1:6]

pgls.out <- list()
for(i in 1:length(pg)){
  mod <- pg[[i]]$model
  sum <- data.frame(summary(mod)$solutions)
  sum$eff.samp <- NULL
  sum <- cbind(rownames(sum), sum)
  rownames(sum) <- NULL
  colnames(sum) <- c('Parameter','mean','lower 95% CI','upper 95% CI','Bayesian p-value')
  sum$mean           <- formatC(sum$mean          , format = "e", digits = 2)
  sum$`lower 95% CI` <- formatC(sum$`lower 95% CI`, format = 'e', digits = 2)
  sum$`upper 95% CI` <- formatC(sum$`upper 95% CI`, format = 'e', digits = 2)
  pgls.out[[i]] <- sum
}
names(pgls.out) <- names(pg)
