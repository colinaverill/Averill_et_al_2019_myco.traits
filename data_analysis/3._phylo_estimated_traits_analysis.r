#Predict traits soley based on phylogeny.
rm(list=ls())
library(picante)
library(data.table)
source('paths.r')
source('functions/phylo_predicted.r')
source('functions/tic_toc.r')

#set output paths.----
out.path <- phylo_estimated_traits.path
full.out.path <- phy_est_models_data.path

#load trait data and phylogeny.----
d <- readRDS(inter_specific_analysis_data.path)
phy <- read.tree(phylogeny_raw.path)
setnames(d,'root_lifespan','root.L')

#subset to observations in phylogeny.----
phy$tip.label <- paste0(toupper(substr(phy$tip.label, 1, 1)), substr(phy$tip.label, 2, nchar(phy$tip.label)))
phy$tip.label <- gsub('_',' ',phy$tip.label)
phy$node.label <- NULL
d <- d[d$Species %in% phy$tip.label,]

#Subset to only species used in final analysis. Must have complete nfix, deciduous, pgf, biome, mat.c, map.c
sub <- d[,c('tpl.Species','mat.c','map.c','pgf','biome3','nfix','deciduous','MYCO_ASSO')]
sub <- sub[complete.cases(sub),]
#Must have at least 1 trait observation.
drop <- d[is.na(d$Ngreen) & is.na(d$Nsenes) & is.na(d$Nroots) & is.na(d$Pgreen) & is.na(d$Psenes) & is.na(d$Proots) & is.na(d$log.LL) & is.na(d$root.L),]
sub <- sub[!(sub$tpl.Species %in% drop$tpl.Species),]
d <- d[d$tpl.Species %in% sub$tpl.Species,]

#Get phylopredicted values.----
tic()
Ngreen <- phylo_predicted(species = d$tpl.Species, trait = log10(d$Ngreen), phy = phy, name = 'Ngreen')
Nsenes <- phylo_predicted(species = d$tpl.Species, trait = log10(d$Nsenes), phy = phy, name = 'Nsenes')
Nroots <- phylo_predicted(species = d$tpl.Species, trait = log10(d$Nroots), phy = phy, name = 'Nroots')
Pgreen <- phylo_predicted(species = d$tpl.Species, trait = log10(d$Pgreen), phy = phy, name = 'Pgreen')
Psenes <- phylo_predicted(species = d$tpl.Species, trait = log10(d$Psenes), phy = phy, name = 'Psenes')
Proots <- phylo_predicted(species = d$tpl.Species, trait = log10(d$Proots), phy = phy, name = 'Proots')
log.LL <- phylo_predicted(species = d$tpl.Species, trait = d$log.LL       , phy = phy, name = 'log.LL')
root.l <- phylo_predicted(species = d$tpl.Species, trait = log10(d$root.L), phy = phy, name = 'root.L')
toc()

#Grab just log10 transformed predicted vs. observed values.----
to_merge <- list(Ngreen, Nsenes, Nroots, Pgreen, Psenes, Proots, log.LL, root.l)
names(to_merge) <- c('Ngreen', 'Nsenes', 'Nroots', 'Pgreen', 'Psenes', 'Proots', 'log.LL', 'root.l')
d <- data.table(d)
d <- d[,.(Species,MYCO_ASSO,Ngreen,Nsenes,Nroots,Pgreen,Psenes,Proots,log.LL,root.L)]
#log10 transform values, save for log.LL
d_keep <- d[,.(Species,MYCO_ASSO,log.LL)]
d_trans <- d[,.(Ngreen,Nsenes,Nroots,Pgreen,Psenes,Proots,root.L)]
d_trans <- apply(d_trans, 2, log10)
d <- data.frame(d_keep,d_trans)
d <- as.data.frame(d)

for(i in 1:length(to_merge)){
  d <- merge(d, to_merge[[i]]$pred, all.x = T)
}

#save trait data with phylogenetic predictions for each observed trait.
saveRDS(d, out.path)
saveRDS(to_merge, full.out.path)
