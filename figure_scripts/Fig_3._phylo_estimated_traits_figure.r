#phylogenetic trait correlation figure.
#clear environment, load packages
rm(list=ls())
source('paths.r')
library(wesanderson)
library(scales)
library(caper)
library(data.table)

#set output path.----
out.path <- phylo_estimated_traits_figure.path
#out.path <- 'test.png'

#load data.----
d <- readRDS(phylo_estimated_traits.path)

#load models, grab lambda values.----
z <- readRDS(phy_est_models_data.path)
lambda <- list()
for(i in 1:length(z)){
  lambda[[i]] <- z[[i]]$lambda
}
lambda <- unlist(lambda)
z <- data.frame(names(z), lambda)
colnames(z)[1] <- 'analysis'


#setup to save.----
png(filename=out.path,width=9.7,height=7.5,units='in',res=300)

#Global plot settings.----
#setup plot
par(mfrow=c(2,3),
    mar = c(1.25,2,1.25,2),
    oma = c(4.5,4.5,3,.2))
#pick colors.
cols <- c('#00acd9','#cfe83c')
r.line.col <- '#ff4ce3'
#assign colors based on mycorrhizal status
d$col <- NA
d$col <- ifelse(d$MYCO_ASSO == 'ECM', cols[2], d$col)
d$col <- ifelse(d$MYCO_ASSO == 'AM' , cols[1], d$col)
#set transparency
trans <- 0.6
#set point size.
p.cex <- 0.8


#Ngreen.----
trait <- 'Ngreen'
label <- 'a. Nitrogen Green'
units <- expression(paste('log(mg N (g tissue)'^'-1',')'))
trait_estimate <- paste0(trait,'_estimate')
tr.form <- formula(paste0(trait,'~',trait_estimate))
mod<- lm(tr.form, data = d)
s.mod <- summary(mod)
#Modified to exclude one very low N green value that throws the whole axis.
plot(tr.form, data = d[d$Ngreen > -1,], cex = p.cex, pch = 16, col = alpha(d$col, trans), xlab = NA, ylab = NA, cex.axis = 1.2, bty ='l')
line.col <- ifelse(s.mod$coefficients[2,4] < 0.05,r.line.col,'gray')
line.type <- ifelse(s.mod$coefficients[2,4] < 0.05, 1,2)
abline(mod, lwd = 2, col = line.col, lty = line.type)
lambda.val <- z[z$analysis == trait,c('lambda')]
lambda.lab <- bquote(lambda == .(format(lambda.val,digits=2)))
mtext(lambda.lab, side = 3, adj = 0.05, line = -7)
mtext(label, side = 3, adj = 0.05, line = -2)
mtext(units, side = 3, adj = 0.05, cex = 0.8, line = -3.5)
mtext(bquote(R^2 == .(round(s.mod$r.squared, 2))), side = 3, adj = 0.05, line = -5.25)

#Nsenes.----
trait <- 'Nsenes'
label <- 'b. Nitrogen Senescent'
units <- expression(paste('log(mg N (g tissue)'^'-1',')'))
mod<- lm(Nsenes ~ Nsenes_estimate, data = d)
s.mod <- summary(mod)
plot(d$Nsenes ~ d$Nsenes_estimate, cex = p.cex, pch = 16, col = alpha(d$col, trans), xlab = NA, ylab = NA, cex.axis = 1.2, bty = 'l')
line.col <- ifelse(s.mod$coefficients[2,4] < 0.05,r.line.col,'gray')
line.type <- ifelse(s.mod$coefficients[2,4] < 0.05, 1,2)
abline(mod, lwd = 2, col = line.col, lty = line.type)
lambda.val <- z[z$analysis == trait,c('lambda')]
lambda.lab <- bquote(lambda == .(format(lambda.val,digits=2)))
mtext(lambda.lab, side = 3, adj = 0.05, line = -7)
mtext(label, side = 3, adj = 0.05, line = -2)
mtext(units, side = 3, adj = 0.05, cex = 0.8, line = -3.5)
mtext(bquote(R^2 == .(round(s.mod$r.squared, 2))), side = 3, adj = 0.05, line = -5.25)

#Nroots.----
trait <- 'Nroots'
label <- 'c. Nitrogen Roots'
units <- expression(paste('log(mg N (g tissue)'^'-1',')'))
trait_estimate <- paste0(trait,'_estimate')
tr.form <- formula(paste0(trait,'~',trait_estimate))
mod<- lm(tr.form, data = d)
s.mod <- summary(mod)
plot(tr.form, data = d, cex = p.cex, pch = 16, col = alpha(d$col, trans), xlab = NA, ylab = NA, cex.axis = 1.2, bty = 'l')
line.col <- ifelse(s.mod$coefficients[2,4] < 0.05,r.line.col,'gray')
line.type <- ifelse(s.mod$coefficients[2,4] < 0.05, 1,2)
abline(mod, lwd = 2, col = line.col, lty = line.type)
lambda.val <- z[z$analysis == trait,c('lambda')]
lambda.lab <- bquote(lambda == .(format(lambda.val,digits=2)))
mtext(lambda.lab, side = 3, adj = 0.05, line = -7)
mtext(label, side = 3, adj = 0.05, line = -2)
mtext(units, side = 3, adj = 0.05, cex = 0.8, line = -3.5)
mtext(bquote(R^2 == .(round(s.mod$r.squared, 2))), side = 3, adj = 0.05, line = -5.25)

#Pgreen.----
trait <- 'Pgreen'
label <- 'd. Phosphorus Green'
units <- expression(paste('log(mg P (g tissue)'^'-1',')'))
trait_estimate <- paste0(trait,'_estimate')
tr.form <- formula(paste0(trait,'~',trait_estimate))
mod<- lm(tr.form, data = d)
s.mod <- summary(mod)
plot(tr.form, data = d, cex = p.cex, pch = 16, col = alpha(d$col, trans), xlab = NA, ylab = NA, cex.axis = 1.2, bty = 'l')
line.col <- ifelse(s.mod$coefficients[2,4] < 0.05,r.line.col,'gray')
line.type <- ifelse(s.mod$coefficients[2,4] < 0.05, 1,2)
abline(mod, lwd = 2, col = line.col, lty = line.type)
lambda.val <- z[z$analysis == trait,c('lambda')]
lambda.lab <- bquote(lambda == .(format(lambda.val,digits=2)))
mtext(lambda.lab, side = 3, adj = 0.05, line = -7)
mtext(label, side = 3, adj = 0.05, line = -2)
mtext(units, side = 3, adj = 0.05, cex = 0.8, line = -3.5)
mtext(bquote(R^2 == .(round(s.mod$r.squared, 2))), side = 3, adj = 0.05, line = -5.25)

#Psenes.----
trait <- 'Pgreen'
label <- 'e. Phosphorus Senescent'
units <- expression(paste('log(mg P (g tissue)'^'-1',')'))
trait_estimate <- paste0(trait,'_estimate')
tr.form <- formula(paste0(trait,'~',trait_estimate))
mod<- lm(tr.form, data = d)
s.mod <- summary(mod)
plot(tr.form, data = d, cex = p.cex, pch = 16, col = alpha(d$col, trans), xlab = NA, ylab = NA, cex.axis = 1.2, bty = 'l')
line.col <- ifelse(s.mod$coefficients[2,4] < 0.05,r.line.col,'gray')
line.type <- ifelse(s.mod$coefficients[2,4] < 0.05, 1,2)
abline(mod, lwd = 2, col = line.col, lty = line.type)
lambda.val <- z[z$analysis == trait,c('lambda')]
lambda.lab <- bquote(lambda == .(format(lambda.val,digits=2)))
mtext(lambda.lab, side = 3, adj = 0.05, line = -7)
mtext(label, side = 3, adj = 0.05, line = -2)
mtext(units, side = 3, adj = 0.05, cex = 0.8, line = -3.5)
mtext(bquote(R^2 == .(round(s.mod$r.squared, 2))), side = 3, adj = 0.05, line = -5.25)

#Proots.----
trait <- 'Proots'
label <- 'f. Phosphorus Roots'
units <- expression(paste('log(mg P (g tissue)'^'-1',')'))
trait_estimate <- paste0(trait,'_estimate')
tr.form <- formula(paste0(trait,'~',trait_estimate))
mod<- lm(tr.form, data = d)
s.mod <- summary(mod)
plot(tr.form, data = d, cex = p.cex, pch = 16, col = alpha(d$col, trans), xlab = NA, ylab = NA, cex.axis = 1.2, bty = 'l')
line.col <- ifelse(s.mod$coefficients[2,4] < 0.05,r.line.col,'gray')
line.type <- ifelse(s.mod$coefficients[2,4] < 0.05, 1,2)
abline(mod, lwd = 2, col = line.col, lty = line.type)
lambda.val <- z[z$analysis == trait,c('lambda')]
lambda.lab <- bquote(lambda == .(format(lambda.val,digits=2)))
mtext(lambda.lab, side = 3, adj = 0.05, line = -7)
mtext(label, side = 3, adj = 0.05, line = -2)
mtext(units, side = 3, adj = 0.05, cex = 0.8, line = -3.5)
mtext(bquote(R^2 == .(round(s.mod$r.squared, 2))), side = 3, adj = 0.05, line = -5.25)

#outer labels.----
mtext('observed trait', side = 2, outer = T, cex = 2, line = 1.5)
mtext('phylogenetic estimated trait', side = 1, outer = T, cex = 2, line = 2.5)
#legend.----
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend('topright',legend = c('arbuscular mycorrhizal','ectomycorrhizal'), 
       pch = 16, col = cols, bty = 'n', cex = 1.6,
       x.intersp = .75, xpd = T, 
       horiz = T)

#end plot.----
dev.off()
