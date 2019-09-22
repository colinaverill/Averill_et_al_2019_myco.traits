#Plot AM-EM trait means across latitudinal zones based on models without variable selection. 
rm(list=ls())
source('paths.r')
source('functions/p_report.r')
library(data.table)
library(caper)
#set output path.----
output.path <- lat_myco_trait_means.path
#output.path <- 'test.png'

#load data.----
d <- readRDS(pgls.glmm_myc.biome3_interaction_no.selection_DECIDUOUS.path)

#subset and order the list.----
list_order <- list('Ngreen','Nsenes','Nroots','Pgreen','Psenes','Proots')
test <- d[names(d) %in% list_order]
names(list_order) <- c('Ngreen','Nsenes','Nroots','Pgreen','Psenes','Proots')
#Get plotting labels.
names(d) <- c('a. Nitrogen Green','b. Nitrogen Senescent','c. Nitrogen Roots',
              'd. Phosphorus Green','e. Phosphorus Senescent','f. Phosphorus Roots')


#setup to save.----
png(filename=output.path,width=9,height=6,units='in',res=300)

#Global plot settings.----
par(mfrow=c(2,3),
    mar = c(1.25,2,1.25,2),
    oma = c(3,4,1,1))
a.cex <- 1.2
trans.box <- c(1,0.8,0.6)
trans.box <- c(0.4,0.3,0.2)
limx <- c(0.8, 2.2)
box.lim <- c(limx[1], limx[1] + 1*(limx[2] - limx[1])/3, limx[1] + 2*(limx[2] - limx[1])/3, limx[2])
inc <- (box.lim[2] - box.lim[1])/2
lat.lab.pos <- c(box.lim[1] + inc,box.lim[2] + inc,box.lim[3] + inc)
x.inter.sp <- 0.08
x <- c(lat.lab.pos[1]-x.inter.sp,lat.lab.pos[1]+x.inter.sp,lat.lab.pos[2]-x.inter.sp,lat.lab.pos[2]+x.inter.sp,lat.lab.pos[3]-x.inter.sp,lat.lab.pos[3]+x.inter.sp)
box.shade.cols <- c('#f9bb41','#9dba32','#2abaac')
myc.col <- c('#f97de2','#9227af')
myc.col <- c('#00acd9','#cfe83c')

#plot loop!----
for(i in 1:6){
  #data unique to model summary output.----
  y <- d[[i]]$mean
  upr <- d[[i]]$upper
  lwr <- d[[i]]$lower
  trait.name <- names(d[i])
  lat.lab <- c('boreal','temperate','tropical')
  bore.N <- d[[i]]$sample_size[c(1,4)]
  temp.N <- d[[i]]$sample_size[c(2,5)]
  trop.N <- d[[i]]$sample_size[c(3,6)]
  N.lab <- c(paste(bore.N, collapse = ', '),
             paste(temp.N, collapse = ', '),
             paste(trop.N, collapse = ', '))
  N.lab <- paste0('N = (',N.lab,')')
  p.val <- d[[i]]$p.values
  p.lab <- unlist(lapply(p.val, p_report))
  
  #adjust margins lower row.
  if(i > 3){par(mar = c(1.25,2,1.25,2))}
  
  #back log transform nutrient data.
  y <- 10^(y)
  upr <- 10^(upr)
  lwr <- 10^(lwr)      
  
  #y limits.
  limy = c(0, max(upr,na.rm = T)*1.05)
  
  #plot box.
  plot(y ~ x, pch = 16, xaxt='n', yaxt = 'n', cex = 0, ylim = limy, xlim = limx, ylab = NA, cex.axis = a.cex, bty='n', xaxs='i', yaxs='i')
  
  #y axis.
  z <- (limy[2] - limy[1])/3
  axis.ticks <- round(c(limy[1], limy[1] + z, limy[1] + 2*z, limy[2]), 1)
  axis(side = 2, at = axis.ticks, las = 2, col = NA, col.ticks = 'black')
  
  #drop rectangles.
  if(i < 4){shade_col <- box.shade.cols[2]}
  if(i > 3){shade_col <- box.shade.cols[3]}
  rect(box.lim[1], limy[1], box.lim[2], limy[2], border = NA, col = adjustcolor(shade_col, trans.box[1]))
  rect(box.lim[2], limy[1], box.lim[3], limy[2], border = NA, col = adjustcolor(shade_col, trans.box[2]))
  rect(box.lim[3], limy[1], box.lim[4], limy[2], border = NA, col = adjustcolor(shade_col, trans.box[3]))
  
  #points and error bars.
  arrows(x, lwr, x, upr, length=0.05, angle=90, code=3)
  #points(y ~ x, pch = 16, col = myc.col, cex = 2)
  points(y ~ x, pch = 21, col = 'black', bg = myc.col, cex = 2)
  
  #latitudinal zone text labels, sample size, significane labels.
  text(lat.lab, x = lat.lab.pos, y = limy[2] * 0.14)
  text(  p.lab, x = lat.lab.pos, y = limy[2] * 0.09)
  text(  N.lab, x = lat.lab.pos, y = limy[2] * 0.03)
  #trait label.
  mtext(trait.name, side = 3, adj = 0.02, line = 0.3)
  
  #unit labels on axes of plots 1 and 4.----
  if(i == 1){mtext(expression(paste('mg N (g tissue)'^'-1')), side = 2, line = 3)}
  if(i == 4){mtext(expression(paste('mg P (g tissue)'^'-1')), side = 2, line = 3)}
}
#Drop legend.----
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('bottom',legend = c('arbuscular mycorrhizal','ectomycorrhizal'), 
       pch = 16, col = myc.col, bty = 'n', cex = 1.6,
       x.intersp = .75, xpd = T, 
       horiz = T)


#end plot.----
dev.off()
