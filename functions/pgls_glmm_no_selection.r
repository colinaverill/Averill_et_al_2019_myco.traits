#' pgls_glmm.r
#' This script is pretty hacky, designed to fit a very particular set of covariates.
#' This script does not perform model selection!
#' Could probably be generalized. 
#'
#' @param pgls.formula   #formula object for fitting models.
#' @param d              #data to fit. Make sure Species is a column in here.
#' @param phy            #your phylogeny.
#'
#' @return               #so much stuff.
#' @export
#'
#' @examples
pgls_glmm <- function(pgls.formula, d, phy){
  
  #compete case the thing.----
  to_keep <- c(all.vars(pgls.formula),'Species','biome3')
  dat <- d[,to_keep]
  dat <- dat[complete.cases(dat),]
  
  #get climate within biome.-----
  temp.mat <- mean(dat[dat$biome3 == 'a_temperate','mat.c'], na.rm = T)
  temp.map <- mean(dat[dat$biome3 == 'a_temperate','map.c'], na.rm = T)
  trop.mat <- mean(dat[dat$biome3 == 'b_tropical' ,'mat.c'], na.rm = T)
  trop.map <- mean(dat[dat$biome3 == 'b_tropical' ,'map.c'], na.rm = T)
  bore.mat <- mean(dat[dat$biome3 == 'c_boreal'   ,'mat.c'], na.rm = T)
  bore.map <- mean(dat[dat$biome3 == 'c_boreal'   ,'map.c'], na.rm = T)
  
  #root-lifespan special case.
  if('root_lifespan' %in% to_keep){
    dat <- dat[dat$biome3 == 'a_temperate',]
  }
  
  #grab sample sizes.----
  dat$biome3 <- substr(dat$biome3,3,nchar(dat$biome3))
  dat$myco_biome <- paste(dat$MYCO_ASSO,dat$biome3)
  sample_size <- table(dat$myco_biome)
  
  #Make phylogenetic random effect (the same as lambda; see Hadfield & Nakagawa 2010).----
  cat('building phylo random effect...\n')
  phy <- force.ultrametric(phy)
  rnd <- inverseA(phy)$Ainv
  #Will Pearse- priors on residual phylogenetic autocorrelation terms.
  priors <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
  
  #fit interactive model.----
  #note- rootlifespan never has an interaction. Hacked to avoid problems with this trait.
  cat('fitting model...\n')
    model <- MCMCglmm(pgls.formula, random=~Species, ginverse=list(Species=rnd), data=dat, prior = priors)
    sum.model <- summary(model)


  #specify design matrix, matched to parameter names.-----
  int <- rep(1,6)
  MYCO_ASSOECM <- c(0,1,0,1,0,1)
  mat.c <- c(bore.mat,bore.mat,temp.mat,temp.mat,trop.mat,trop.mat)
  map.c <- c(bore.map,bore.map,temp.map,temp.map,trop.map,trop.map)
  biome_bore <- c(1,1,0,0,0,0)
  biome_trop <- c(0,0,0,0,1,1)
  MYCO_ASSOECM_biome_trop <- c(0,0,0,0,0,1)
  MYCO_ASSOECM_biome_bore <- c(0,1,0,0,0,0)
  design.matrix <- cbind(int,MYCO_ASSOECM,mat.c,map.c,biome_bore,biome_trop,MYCO_ASSOECM_biome_trop,MYCO_ASSOECM_biome_bore)
  colnames(design.matrix)[1] <- '(Intercept)'
  colnames(design.matrix)[7] <- 'MYCO_ASSOECM:biome_trop'
  colnames(design.matrix)[8] <- 'MYCO_ASSOECM:biome_bore'
  
  #get posteriors of mycorrhizal types by biome.----
  mcmc <- model$Sol
  group.list <- list()
  for(i in 1:10000){
    par <- mcmc[sample(nrow(mcmc),1),]
    this.matrix <- design.matrix[,colnames(design.matrix) %in% names(par)]
    par <- par[names(par) %in% colnames(this.matrix)]
    this.matrix <- this.matrix[,order(match(colnames(this.matrix), names(par)))]
    output <- as.vector(this.matrix %*% par)
    names(output) <- c('am_bore','em_bore','am_temp','em_temp','am_trop','em_trop')
    group.list[[i]] <- output
  }
  group.list <- do.call(rbind,group.list)
  
  #get mean and 95% CI of mycorrhizal contrasts within biome, as well as p-values.----
  #mean
  mu <- colMeans(group.list)
  #contrast of within biome mycorrhizal posteriors.
  check.temp <- group.list[,4] - group.list[,3]
  check.trop <- group.list[,6] - group.list[,5]
  check.bore <- group.list[,2] - group.list[,1]
  #Get fraction of parameter draws that are above zero.
  p.temp <- sum(check.temp < 0) / length(check.temp)
  p.trop <- sum(check.trop < 0) / length(check.trop)
  p.bore <- sum(check.bore < 0) / length(check.bore)
  #if p > 0.5, means more parameter density below zero. Subtract p-value from 1.
  if(p.temp > 0.5){p.temp = 1-p.temp}
  if(p.trop > 0.5){p.trop = 1-p.trop}
  if(p.bore > 0.5){p.bore = 1-p.bore}
  #two tailed test, multiply computed p-values by two.
  p.temp <- p.temp * 2
  p.trop <- p.trop * 2
  p.bore <- p.bore * 2
  p.values <- c(p.bore,p.temp,p.trop)
  #get posterior contrast 95% CI.
  check.temp <- quantile(check.temp, prob = c(0.025, 0.975))
  check.trop <- quantile(check.trop, prob = c(0.025, 0.975))
  check.bore <- quantile(check.bore, prob = c(0.025, 0.975))
  #calculate standard error. CI is mu +/- 2 SE. 
  se.temp <- abs(check.temp[1] - check.temp[2]) / 4
  se.trop <- abs(check.trop[1] - check.trop[2]) / 4
  se.bore <- abs(check.bore[1] - check.bore[2]) / 4
  #assign significance.
  if(check.temp[1] * check.temp[2] > 0){signif.temp = '*'}
  if(check.trop[1] * check.trop[2] > 0){signif.trop = '*'}
  if(check.bore[1] * check.bore[2] > 0){signif.bore = '*'}
  if(exists('signif.temp') == F){signif.temp = ''}
  if(exists('signif.trop') == F){signif.trop = ''}
  if(exists('signif.bore') == F){signif.bore = ''}
  signif <- c(signif.bore, signif.temp, signif.trop)
  se <- c(se.bore, se.bore, se.temp, se.temp, se.trop, se.trop)
  lo95 <- mu - se
  hi95 <- mu + se
  
  #return output.
  output <- list(model, dat, group.list,mu,se,lo95,hi95,signif,p.values,pgls.formula,sample_size)
  names(output) <- c('model','data','posteriors','mean','se','lower','upper','signif','p.values','final_model_formula','sample_size')
  return(output)  
}
