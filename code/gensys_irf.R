library(tsibble)
library(scales)

gensys_irf_ts <- function(gsys, shocks, nirf, varnames, shocknames) {
  
  nshocks <- length(shocks)
  nvar <- length(varnames)
  ndf <- nshocks*nvar*nirf
    
  irf.df <- tibble(`Shock`=NA, `Variable`=NA, `Period`=NA, `Impact`=NA)
  
  M <- gsys$impact_sol
  G <- gsys$G_sol
  #C <- gsys$cons_sol
  
  irfidx <-1
  
  for(sidx in 1:length(shocks)) {
    shock_vector <- matrix(nrow=length(shocks), ncol=1, data=0)
    shock_vector[sidx] <- shocks[sidx]
    
    impact0 <- M %*% shock_vector
    sim <- 1
    for(v in 1:nvar) {
      irf.df[irfidx,"Shock"] <- shocknames[sidx]
      irf.df[irfidx, "Variable"] <- varnames[v]
      irf.df[irfidx, "Period"] <- sim
      irf.df[irfidx, "Impact"] <- impact0[v,1] 
      irfidx <- irfidx + 1
    }
    
    for(sim in 2:nirf) {
      impact1 <- G %*% impact0
      for(v in 1:nvar) {
        irf.df[irfidx,"Shock"] <- shocknames[sidx]
        irf.df[irfidx, "Variable"] <- varnames[v]
        irf.df[irfidx, "Period"] <- sim
        irf.df[irfidx, "Impact"] <- impact1[v,1] 
        irfidx <- irfidx + 1
      }
      impact0 <- impact1
    }
  }
  
  irf.ts <- as_tsibble(irf.df, key=c(Shock, Variable), index=Period)
  irf.ts %>%
    mutate(Shock = factor(Shock, levels=shocknames, ordered=TRUE)) %>%
    mutate(Variable = factor(Variable, levels=varnames, ordered=TRUE)) ->
    irf.ts

  return(irf.ts)
}

plot_irf <- function(irf.ts) {
  ggplot(irf.ts, aes(x=Period, y=Impact)) +
    geom_line(size=1.2, color="dodgerblue4") +
    facet_grid(Variable~Shock, scales="free_y") +
    theme_bw() +
    scale_x_continuous(breaks=pretty_breaks())-> 
    plt
  
  return(plt)
}
