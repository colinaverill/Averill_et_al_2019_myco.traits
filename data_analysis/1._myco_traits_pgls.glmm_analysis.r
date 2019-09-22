#Using MCMCglmm to fit PGLS models with deciduousness as a covariate.
#zero model selection performed, just straight up sending model structures based on apriori hypotheses.
rm(list=ls())
source('paths.r')
source('functions/pgls_glmm_no_selection.r')
source('functions/tic_toc.r')
library(data.table)
library(phytools)
library(MCMCglmm)
library(caper)

#set output path.----
output.path <- pgls.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.path

#load data.----
d <- readRDS(inter_specific_analysis_data.path)
phy <- read.tree(phylogeny_raw.path) #'colin_2018-12--2.tre'
d$biome_trop <- ifelse(d$biome3 == 'b_tropical',1, 0)
d$biome_bore <- ifelse(d$biome3 == 'c_boreal',1,0)

#Some data prep.----
phy$tip.label <- paste0(toupper(substr(phy$tip.label, 1, 1)), substr(phy$tip.label, 2, nchar(phy$tip.label)))
phy$tip.label <- gsub('_',' ',phy$tip.label)
phy$node.label <- NULL
d <- d[d$Species %in% phy$tip.label,]

#Specify model structure (predictors)----
#unfortunately, the spaces are important.
predictors <- 'MYCO_ASSO + MYCO_ASSO:biome_bore + MYCO_ASSO:biome_trop + biome_trop + biome_bore + nfix + pgf + deciduous + mat.c + map.c'

#specify traits.----
traits <- c('Ngreen','Nsenes','Nroots','Pgreen','Psenes','Proots')

#fit your models.----
all.fit <- list()
tic()
for(i in 1:length(traits)) {
  form <- paste0('log10(',traits[i],') ~ ',predictors)
  form <- as.formula(form)
  fit <- pgls_glmm(form, d, phy)
  all.fit[[i]] <- fit
  cat(traits[i],'model fitted.',i,'of',length(traits),'complete. ');toc()
}
names(all.fit) <- traits

#save output.----
saveRDS(all.fit, output.path)
