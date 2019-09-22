#report p values for figures.
p_report <- function(x){
  if(x >= 0.10)          {p = paste0('N.S.')}
  if(x >= 0.05 & x < 0.1){p = paste0('p = ',round(x,2))}
  if(x < 0.05 & x > 0.01){p = paste0('p = ',round(x,2))}
  if(x < 0.01)           {p = 'p < 0.01'}
  return(p)
}
