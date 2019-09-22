#' phylo_predicted: this function takes a phylogeny and set of known trait values. 
#' then estimates each trait value from the phylogeny by sequentially leaving one observation out and making a phylogenetic prediction using phyEstimate from the picante package.
#' this function depends on the picante package.
#'
#' @param species vector of species names.
#' @param trait vector of corresponding trait values (can contain NAs)
#' @param phy phylogeny. Can contain species without trait observations. species names must match phylogeny tip labels.
#' @param name optional. Name of traits to return a data frame with named parameter columns.
#'
#' @return returns a dataframe with phylogenetic estimates of trait values for each species, a se estimate, and a species name column.
#' @export
#'
#' @examples
#' library(picante)
#' phy <- read.tree('/home/caverill/FIA_trait_analyses/raw_data/SuperTree_2017_09_26')
#' traits <- read.csv('/home/caverill/FIA_trait_analyses/raw_data/species_traitsCheckedPhylo_2017_08_16.csv')
#' test <- phylo_predicted(traits$species,traits$NitrogenGreen, phy, name = 'Ngreen')
#'
phylo_predicted <- function(species, trait, phy, name=NA){
  source('functions/pic_pro.r')
  library(MCMCglmm)
  library(phytools)
  
  #grab complete observations of traits, assign names.
  dat <- data.frame(species, trait)
  dat <- dat[complete.cases(dat),]
  dat$species <- as.character(dat$species)
        trait  <- dat$trait
  names(trait) <- dat$species
  
  #make sure phylogeny is happy.
  phy$node.label <- NULL
  phy$tip.label  <- paste0(toupper(substr(phy$tip.label, 1, 1)), substr(phy$tip.label, 2, nchar(phy$tip.label)))
  phy$tip.label <-  gsub('_',' ',phy$tip.label) 
  
  #prune the phylogeny to match.
  to.drop <- phy$tip.label[!(phy$tip.label %in% dat$species)]
  n.phy <- ape::drop.tip(phy,to.drop)
  n.phy <- multi2di(n.phy) #deal with polytomies.
  
  #match the trait file again.
  trait <- trait[names(trait) %in% n.phy$tip.label]
  
  #run phyEstimate, leaving one observation out at a time.
  trait_est <- list()
  for(i in seq_along(trait)){
    test.trait <- trait[-i]
    trait_est[[i]] <- picante::phyEstimate(n.phy, test.trait)
  }
  #collapse list, return as data.frame
  trait_est <- do.call('rbind',trait_est)
  trait_est$Species <- rownames(trait_est)
  #name the trait estimate in output.
  if(!is.na(name)){colnames(trait_est) <- c(paste0(name,'_estimate'), paste0(name,'_se'),'Species')}
  
  #Get a lambda value for the intercept only model.
  phy <- force.ultrametric(phy)
  rnd <- inverseA(phy)$Ainv
  #priors on residual phylogenetic autocorrelation terms.
  priors <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
  colnames(dat)[2] <- 'trait.1' #MCMCglmm will not accept 'trait' as a variable, cause its special.
  mod <- MCMCglmm(trait.1 ~ 1, random=~species, ginverse=list(species=rnd), data=dat, prior = priors)
  lambda <- mod$VCV[,'species']/(mod$VCV[,'species']+mod$VCV[,'units'])
  lambda <- mean(lambda)

  #create list of output
  output <- list(trait_est,mod,lambda)
  names(output) <- c('pred','model','lambda')

  #return output.
  return(output)
}