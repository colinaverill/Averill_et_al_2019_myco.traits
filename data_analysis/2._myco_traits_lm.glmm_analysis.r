#Using MCMCglmm to fit lm models with same structure as selected pgls models.
rm(list=ls())
source('paths.r')
source('functions/lm_glmm.r')
source('functions/tic_toc.r')
library(data.table)
library(phytools)
library(MCMCglmm)
library(caper)

#set output path.----
output.path <- lm.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.path

#load data and pgls results.----
pgls.fits <- readRDS(pgls.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.path)
ref.dat <- readRDS(inter_specific_analysis_data.path)
ref.dat$biome_trop <- ifelse(ref.dat$biome3 == 'b_tropical',1, 0)
ref.dat$biome_bore <- ifelse(ref.dat$biome3 == 'c_boreal'  ,1, 0)


#fit your models.----
all.fit <- list()
tic()
for(i in 1:length(pgls.fits)) {
  #fit.
  all.fit[[i]] <- lm_glmm(pgls.fits[[i]]$final_model_formula, ref.dat)
  #report.
  cat(names(pgls.fits)[i],'model fitted.',i,'of',length(pgls.fits),'complete. ');toc()
}
names(all.fit) <- names(pgls.fits)

#save output.----
saveRDS(all.fit, output.path)
