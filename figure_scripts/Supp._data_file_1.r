#Generating supplementary data file 1. Species list with myco associations, plant growth form and N fixation
rm(list=ls())
source('paths.r')
library(caper)

#set output path.----
output.path <- Supplementary_Data_File_1.path

#Data subsetting.----
#Filter based on interspecific observations actually used in the analysis.
inter <- readRDS(inter_specific_analysis_data.path)
phy <- read.tree(phylogeny_raw.path) #'colin_2018-12--2.tre'

#Some data manipulation so I can match intra-specific observations to speciesa included in interspecific analysis.
phy$tip.label <- paste0(toupper(substr(phy$tip.label, 1, 1)), substr(phy$tip.label, 2, nchar(phy$tip.label)))
phy$tip.label <- gsub('_',' ',phy$tip.label)
phy$node.label <- NULL
inter <- inter[inter$Species %in% phy$tip.label,]
#Must have at least 1 trait observation.
drop <- inter[is.na(inter$Ngreen) & is.na(inter$Nsenes) & is.na(inter$Nroots) & is.na(inter$Pgreen) & is.na(inter$Psenes) & is.na(inter$Proots) & is.na(inter$log.LL) & is.na(inter$root_lifespan),]
inter <- inter[!(inter$tpl.Species %in% drop$tpl.Species),]
inter <- inter[,c('tpl.Species','biome3','MYCO_ASSO','nfix','pgf','mat.c','map.c','deciduous','myco_doi')]
inter <- inter[complete.cases(inter),]
inter$mat.c <- NULL
inter$map.c <- NULL

#clean up the names.----
out <- inter
colnames(out) <- c('species','latitudinal_zone','mycorrhizal_association','nitrogen_fixation_status','plant_growth_form','leaf_habit','mycorrhizal_reference')
out$latitudinal_zone <- substr(out$latitudinal_zone,3,nchar(out$latitudinal_zone))
out$plant_growth_form <- ifelse(out$plant_growth_form == 'angio','angiosperm',out$plant_growth_form)
out$plant_growth_form <- ifelse(out$plant_growth_form == 'gymno','gymnosperm',out$plant_growth_form)
out$leaf_habit        <- ifelse(out$leaf_habit == 0, 'evergreen','deciduous')
#Save output.----
write.csv(out, output.path)
