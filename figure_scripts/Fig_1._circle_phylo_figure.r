#Phylogeny figure.
rm(list=ls())
source('paths.r')
library(caper)
library(scales)
library(data.table)
library(willeerd)

#load trait data and phylogeny.----
d <- readRDS(inter_specific_analysis_data.path)
phy <- read.tree(phylogeny_raw.path)
phy$tip.label  <- paste0(toupper(substr(phy$tip.label, 1, 1)), substr(phy$tip.label, 2, nchar(phy$tip.label)))
phy$tip.label <- gsub("_", " ", phy$tip.label)
phy <- drop.tip(phy, setdiff(phy$tip.label, d$Species))

#only include species in phylogeny.
d <- d[d$Species %in% phy$tip.label,]

myc.assoc <- with(d, setNames(MYCO_ASSO, Species))
myc.assoc <- as.character(myc.assoc[phy$tip.label])
pgf.assoc <- with(d, setNames(      pgf, Species))
pgf.assoc <- as.character(pgf.assoc[phy$tip.label])

#pick colors
cols1 <- c('#00acd9','#cfe83c') #AM, ECM colors
cols2 <- c('#ffd821','#ff8f00') #gymno angio colors
cols2 <- c('purple','#ffd821') #gymno angio colors
#set transparency
trans <- 0.5


#Begin plot commands.----
png(filename=phylogeny_figure.path,width=8,height=6,units='in',res=300)

par(mar = c(2,13,2,1), xpd = NA)
cartoon.plot(phy, type="fan", show.tip.label=FALSE, auto.polies=TRUE, clade.col="grey40")
tipring(tips=which(myc.assoc=="AM" ), col=alpha(cols1[1], trans), lwd=8, radial.adj=1.1)
tipring(tips=which(myc.assoc=="ECM"), col=alpha(cols1[2], trans), lwd=8, radial.adj=1.1)
tipring(tips=which(pgf.assoc=='gymno'), col = cols2[1], lwd = 6, radial.adj = 1.05, pch=15)
tipring(tips=which(pgf.assoc=='angio'), col = cols2[2], lwd = 6, radial.adj = 1.05, pch=15)
#mycorrhizal legend.
mtext(expression(paste(bold('mycorrhizal association'))), line = -5.1, adj = -.725)
legend(-800,275, legend = c('arbuscular \nmycorrhizal','ectomycorrhizal'), 
       pch = 16, col = c(cols1[1],cols1[2]), 
       bty = 'n', x.intersp = 0.7, y.intersp = 1.5, cex = 1)
mtext(expression(paste(bold('plant growth form'))), line = -13.5, adj = -.625)
legend(-800,-000, legend = c('gymnosperm','angiosperm'), 
       pch = 16, col = c(cols2[1],cols2[2]), 
       bty = 'n', x.intersp = 0.7, y.intersp = 1.5, cex = 1)

#end plot.
dev.off()
