#Plotting lm vs. PGLS results of tissue chemistry traits.
#clear envrionment, load packages.
rm(list=ls())
library(data.table)
source('paths.r')
source('functions/p_report.r')

#set output path.----
output.path <- lm_pgls_effects_figure.path
#output.path <- 'test.png'

#load summary data for figure.----
pg <- readRDS(pgls.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.path)
lm <- readRDS(  lm.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.path)

#subset and order the list.----
list_order     <- list('Ngreen','Nsenes','Nroots','Pgreen','Psenes','Proots')
names(list_order) <- c('Ngreen','Nsenes','Nroots','Pgreen','Psenes','Proots')
pg <- pg[names(list_order)]
lm <- lm[names(list_order)]
#Get plotting labels.
to_name <- c('Nitrogen Green','Nitrogen Senescent','Nitrogen Roots',
             'Phosphorus Green','Phosphorus Senescent','Phosphorus Roots')
to_name <- c('a. Nitrogen Green','b. Nitrogen Senescent','c. Nitrogen Roots',
             'd. Phosphorus Green','e. Phosphorus Senescent','f. Phosphorus Roots')
names(pg) <- to_name
names(lm) <- to_name

#setup to save.----
png(filename=output.path,width=9,height=7,units='in',res=300)

#global plot settings.-----
#setup plot panels and margins.
par(mfrow=c(2,3),
    mar = c(1,2,1,2),
    oma = c(5,5,3,1))
#get x positions
x <- c(0.9,1.05,1.4,1.55)
limx <- c(0.8,1.7)
a.cex <- 1.2
#pick colors.
cols <- c('#00acd9','#cfe83c')

#begin plot loop.----
for(i in 1:length(pg)){
  mean <- 10^c(lm[[i]]$mean[3:4]  , pg[[i]]$mean[3:4]  )
  lo95 <- 10^c(lm[[i]]$lower[3:4] , pg[[i]]$lower[3:4] )
  hi95 <- 10^c(lm[[i]]$upper[3:4] , pg[[i]]$upper[3:4] )
  p <-    c(lm[[i]]$p.values[2], pg[[i]]$p.values[2])
  limy <- c(0,max(mean)*1.5)
  plot(mean ~ x, pch = 16, xaxt='n', cex = 0, ylim = limy, xlim = limx, ylab = NA, col = cols, cex.axis = a.cex, bty = 'l')
  arrows(x, lo95, x, hi95, length=0.05, angle=90, code=3)
  points(mean ~ x, pch = 16, col = cols, cex = 2)
  abline(v=((limx[2] - limx[1])/2) + limx[1], lwd = 2, lty = 2, col = 'light gray')
  mtext(names(pg)[i], side = 3, adj = 0.05, line = -2)
  mtext(p_report(p[1]), side = 1, line = -1.5, adj = 0.15)
  mtext(p_report(p[2]), side = 1, line = -1.5, adj = 0.8)
  #units
  if(i == 1){mtext(expression(paste('mg N (g tissue)'^'-1')), side = 2, cex = 1.3, line = 2.5)}
  if(i == 4){mtext(expression(paste('mg P (g tissue)'^'-1')), side = 2, cex = 1.3, line = 2.5)}
}

#outer labels.----
mtext('without', side = 1, outer = T, cex = 1.2, line = .5, adj = 0.07)
mtext('correction', side = 1, outer = T, cex = 1.2, line = 2, adj = 0.05)
mtext('phylo-', side = 1, outer = T, cex = 1.2, line =.5, adj = 0.21)
mtext('corrected', side = 1, outer = T, cex = 1.2, line = 2, adj = 0.20)

mtext('without', side = 1, outer = T, cex = 1.2, line = .5, adj = 0.43)
mtext('correction', side = 1, outer = T, cex = 1.2, line = 2, adj = 0.43)
mtext('phylo-', side = 1, outer = T, cex = 1.2, line =.5, adj = 0.57)
mtext('corrected', side = 1, outer = T, cex = 1.2, line = 2, adj = 0.57)

mtext('without', side = 1, outer = T, cex = 1.2, line = .5, adj = 0.79)
mtext('correction', side = 1, outer = T, cex = 1.2, line = 2, adj = 0.80)
mtext('phylo-', side = 1, outer = T, cex = 1.2, line =.5, adj = 0.94)
mtext('corrected', side = 1, outer = T, cex = 1.2, line = 2, adj = 0.95)

#Drop legend.----
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('top',legend = c('arbuscular mycorrhizal','ectomycorrhizal'), 
       pch = 16, col = cols, bty = 'n', cex = 1.6,
       x.intersp = .75, xpd = T, 
       horiz = T)


#end plot.----
dev.off()