#Plotting global distribution of trait observations used in this analysis.
rm(list=ls())
source('paths.r')
library(ggplot2)
library(ggalt)
library(data.table)

#set output path.----
output.path <- sample_map_figure.path
ouput.path <- 'map.png'

#grab and filter interspecific data for subsetting.----
d <- readRDS(inter_specific_analysis_data.path)
phy <- read.tree(phylogeny_raw.path) #'colin_2018-12--2.tre'
d$biome_trop <- ifelse(d$biome3 == 'b_tropical',1, 0)
d$biome_bore <- ifelse(d$biome3 ==   'c_boreal',1, 0)

#Some data prep.----
phy$tip.label <- paste0(toupper(substr(phy$tip.label, 1, 1)), substr(phy$tip.label, 2, nchar(phy$tip.label)))
phy$tip.label <- gsub('_',' ',phy$tip.label)
phy$node.label <- NULL
d <- d[d$Species %in% phy$tip.label,]

#Must have at least 1 trait observation.
drop <- d[is.na(d$Ngreen) & is.na(d$Nsenes) & is.na(d$Nroots) & is.na(d$Pgreen) & is.na(d$Psenes) & is.na(d$Proots) & is.na(d$log.LL) & is.na(d$root_lifespan),]
d <- d[!(d$tpl.Species %in% drop$tpl.Species),]
preds <- c('tpl.Species','MYCO_ASSO','nfix','pgf','deciduous','mat.c','map.c','biome_bore','biome_trop')
keep <- d[,preds]
keep <- keep[complete.cases(keep),]
d <- d[d$tpl.Species %in% keep$tpl.Species,]
inter <- d

#load intra-specific data.----
d <- data.table(readRDS(intra_specific_analysis_data.path))
d <- d[d$tpl.Species %in% inter$tpl.Species,]
d <- d[,.(latitude, longitude)]
d <- d[complete.cases(d),]
lon <- d$longitude
lat <- d$latitude

#set output spec.----
png(filename=output.path,width=8,height=4.8,units='in',res=300)

world <- map_data('world')
map <- ggplot() + geom_cartogram(data = world, map = world, 
                                 aes(x=long, y = lat, group = group, map_id=region))
map <- map + coord_proj("+proj=wintri", ylim = c(-55,90))
map <- map + geom_point(aes(x = lon, y = lat), color = "yellow"    , size = .2)
map <- map + theme(axis.line=element_blank(),
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   legend.position="none",
                   panel.border=element_blank()
                  )
map

#end plot.----
dev.off()