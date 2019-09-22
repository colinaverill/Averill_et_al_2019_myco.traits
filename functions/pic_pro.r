#this function takes a dataframe, a phylogeny, and some traits.
#trait.1 is the dependent varuiable.
#trait.2 - trait.N are predictors to be included in PGLS model.
#This function prunes the phylogeny, builds a comparative data object, and fits a PGLS model.

pic_pro <- function(y,x, phylogeny, trait.data, intercept = T, int.only = F, log = F, interaction = NA, lambda = 'ML'){
  #make sure you aren't dealing with data.table.
  trait.data <- as.data.frame(trait.data)

  #subset data to only include species with complete data
  to_keep <- c(y,x)
  to_keep <- to_keep[!is.na(to_keep)]
  d.sub <- trait.data[,c(to_keep,'Species')]
  d.sub <- d.sub[complete.cases(d.sub),]
  
  #drop levels.
  d.sub <- droplevels.data.frame(d.sub)
  
  #prune the phylogeny to match.
  #check for underscores, make sure everything the same.
  if(length(grep('_',phylogeny$tip.label)) > 0){
    d.sub$Species <- gsub(' ','_',d.sub$Species)
    cat('Underscores present in phylogeny tip labels. Modifying species labels to match.\n')
  }
  to.drop <- phylogeny$tip.label[!(phylogeny$tip.label %in% d.sub$Species)]
  n.phy <- drop.tip(phylogeny,to.drop)
  
  #make a comparative data object for PGLS in caper.
  c.data <- comparative.data(n.phy, d.sub, Species)
  
  #Cook up formula for PGLS
  if(int.only == T){
    linear <- '1'
  }
  if(int.only == F){
    linear <- paste0(paste0(to_keep[2:length(to_keep)], sep = '+'), collapse = '')
    linear <- substr(linear,1,nchar(linear)-1) #remove trailing '+'
    if(!is.na(interaction)){
      insert <- paste(interaction, collapse = '+') 
      linear <- paste0(linear, '+', insert)
    }
  }
  formula <- paste((y),'~')
  if(log == T){
    formula <- paste0('log10(',y,') ~ ')
    if(y == 'log.LL'){
      formula <- paste((y),'~')
    }
  }
  formula <- paste(formula,linear)
  if(intercept == F){formula <- paste(formula,' - 1')}
  formula <- as.formula(formula)
  
  #run PGLS, return output
  model <- try(pgls(formula, data = c.data, lambda = lambda))
  #Sometimes the above fails. If that's the case just set lambda to 1.
  if(class(model) == 'try-error'){
    model <- pgls(formula, data = c.data, lambda = 1)
  }
  return(model)
} #end function.
