library(tidyverse)
library(BMR)
library(QZ)

source("gensys_irf.R")
#devtools::install_github("kthohr/BMR")

# Put New Keynesian into Sims's GenSys form: Gamma_O y_t = C + Gamma_1 y_{t-1} + \Psi \epsilon_t + \Pi \eta_t
# where y_t = E_{t-1} y_{t-1} + \eta_t

# Indices into matrix y_t
i_var_x <- 1        # Output gap x_t
i_var_pi <- 2       # Inflation rate \pi_t
i_var_r <- 3        # Federal funds rate r_t
i_var_xe <- 4       # Expectation E_t x_{t+1}
i_var_pie <- 5      # Expectation E_t \pi_{t+1}
i_var_shock_x <- 6  # Demand shock
i_var_shock_pi <- 7 # Inflation shock
i_var_shock_r <- 8  # Monetary policy shock
i_var_piA <- 9      # Inflation target, i.e. some average
i_var_piB <- 10     # Backward-looking target
i_var_piF <- 11     # Forward-looking target
i_var_piFe <- 12    # Expectation of Forward-looking target
i_var_xe_alt <- 13  # Non-rational expectation for future output gap
i_var_pie_alt <- 14 # Non-rational expectation for future inflation
nvar <- 14

varnames <- rep(NA, nvar)
varnames[i_var_x] <- "Output Gap"
varnames[i_var_pi] <- "Inflation"
varnames[i_var_r] <- "Federal Funds Rate"
varnames[i_var_xe] <- "Rational Expected Output Gap"
varnames[i_var_pie] <- "Rational Expected Inflation Rate"
varnames[i_var_shock_x] <- "Demand Shock"
varnames[i_var_shock_pi] <- "Cost Shock"
varnames[i_var_shock_r] <- "Monetary Policy Rate"
varnames[i_var_piA] <- "Inflation Target"
varnames[i_var_piB] <- "Backward-Looking Average Inflation"
varnames[i_var_piF] <- "Forward-Looking Average Inflation"
varnames[i_var_piFe] <- "Rational Expected Forward-Looking Average Inflation"
varnames[i_var_xe_alt] <- "Non-Rational Expectation for Output Gap" 
varnames[i_var_pie_alt] <- "Non-Rational Expectation for Inflation" 

# Shock indices
i_shock_x <- 1
i_shock_pi <- 2
i_shock_r <- 3
nshocks <- 3

shocknames <- rep(NA, nshocks)
shocknames[i_shock_x] <- "Demand Shock"
shocknames[i_shock_pi] <- "Cost Shock"
shocknames[i_shock_r] <- "Monetary Policy Shock"

# Expectation error indices
i_exp_xe <- 1
i_exp_pie <- 2
i_exp_piFe <- 3
nexp <- 3

# Equation indices
eq_IS <- 1
eq_Phillips <- 2
eq_Taylor <- 3
eq_xe <- 4
eq_pie <- 5
eq_shock_x <- 6
eq_shock_pi <- 7
eq_shock_r <- 8
eq_piA <- 9
eq_piB <- 10
eq_piF <- 11
eq_piFe <- 12
eq_xe_alt <- 13
eq_pie_alt <- 14
neq <- 14

nksys <- function(gamma, deltaB, deltaF, lambda)
{
  # Setup matrices
  Gamma0 <- matrix(nrow=neq, ncol=nvar, data=0)
  C <- matrix(nrow=neq, ncol=1, data=0)
  Gamma1 <- matrix(nrow=neq, ncol=nvar, data=0)
  Psi <- matrix(nrow=neq, ncol=nshocks, data=0)
  Pi <- matrix(nrow=neq, ncol=nexp, data=0)
  Sigma <- matrix(nrow=nshocks, ncol=nshocks, data=0)
  
  
  # IS equation
  # x_t = x_{t+1}^e - 1/sigma \left[(r_t - r*) - (\pi_{t+1}^e - \pi*)\right] + \epsilon_t^x
  Gamma0[eq_IS, i_var_x] <- 1.0
  Gamma0[eq_IS, i_var_r] <- 1.0/sigma
  Gamma0[eq_IS, i_var_xe_alt] <- -1.0 
  Gamma0[eq_IS, i_var_pie_alt] <- -1.0/sigma
  C[eq_IS, 1] <- 1.0/sigma*(rstar-pistar)
  Gamma0[eq_IS, i_var_shock_x] <- -1.0
  
  # Phillips curve
  # (\pi_t - \pi*) = \beta (\pi_{t+1}^e - \pi*) + \kappa x_t + \epsilon_t^\pi
  Gamma0[eq_Phillips, i_var_x] <- -1.0*kappa
  Gamma0[eq_Phillips, i_var_pi] <- 1.0
  Gamma0[eq_Phillips, i_var_pie_alt] <- -1.0*beta
  C[eq_Phillips, 1] <- (1.0-beta)*pistar
  Gamma0[eq_Phillips, i_var_shock_pi] <- -1.0
  
  # Taylor Rule
  # r_t - r* = \rho_r (r_{t-1} - r*) + (1-\rho_r) \left[ \psi_\pi (\pi_t^A - \pi^*) + \psi_x x_t \right] + \epsilon_t^r
  Gamma0[eq_Taylor, i_var_x] <- -1.0*(1.0-rho_r)*psi_x
  Gamma0[eq_Taylor, i_var_piA] <- -1.0*(1.0-rho_r)*psi_pi
  Gamma0[eq_Taylor, i_var_r] <- 1.0
  C[eq_Taylor, 1] <- (1.0-rho_r)*(rstar - psi_pi*pistar)
  Gamma1[eq_Taylor, i_var_r] <- rho_r
  Gamma0[eq_Taylor, i_var_shock_r] <- -1.0
  
  # Average inflation
  # \pi_t^A = \gamma \pi_t^B + (1 - \gamma) \pi_t^F
  Gamma0[eq_piA, i_var_piA] = 1.0
  Gamma0[eq_piA, i_var_piB] = -1.0*gamma
  Gamma0[eq_piA, i_var_piF] = -1.0*(1.0 - gamma)
  
  # Backward-looking average inflation
  # \pi_t^B = (1 - \delta_B) \pi_{t-1}^B + \delta_B \pi_t
  Gamma0[eq_piB, i_var_piB] <- 1.0
  Gamma1[eq_piB, i_var_piB] <- 1.0 - deltaB
  Gamma0[eq_piB, i_var_pi] <- -1.0*deltaB
  
  #Forward-looking average inflation
  # \pi_t^F = (1 - \delta_F) E_t \pi_{t+1}^F + \delta_F E_t \pi_{t+1}
  Gamma0[eq_piF, i_var_piF] <- 1.0
  Gamma0[eq_piF, i_var_piFe] <- -1.0*(1.0-deltaF)
  Gamma0[eq_piF, i_var_pie] <- -1.0*deltaF
  
  
  # Output Gap Expectation
  # x_t = E_{t-1} x_t + \eta^x,t 
  Gamma0[eq_xe, i_var_x] <- 1.0
  Gamma1[eq_xe, i_var_xe] <- 1.0
  Pi[eq_xe, i_exp_xe] <- 1.0
  
  # Inflation Expectation
  Gamma0[eq_pie, i_var_pi] <- 1.0
  Gamma1[eq_pie, i_var_pie] <- 1.0
  Pi[eq_pie, i_exp_pie] <- 1.0
  
  # Forward Average Inflation Expectation
  Gamma0[eq_piFe, i_var_piF] <- 1.0
  Gamma1[eq_piFe, i_var_piFe] <- 1.0
  Pi[eq_piFe, i_exp_piFe] <- 1.0
  
  # Demand Shock
  Gamma0[eq_shock_x, i_var_shock_x] <- 1.0
  Gamma1[eq_shock_x, i_var_shock_x] <- rho_x
  Psi[eq_shock_x, i_shock_x] <- 1.0
  
  # Cost Shock
  Gamma0[eq_shock_pi, i_var_shock_pi] <- 1.0
  Gamma1[eq_shock_pi, i_var_shock_pi] <- rho_pi
  Psi[eq_shock_pi, i_shock_pi] <- 1.0
  
  # Monetary Policy Shock
  Gamma0[eq_shock_r, i_var_shock_r] <- 1.0
  Psi[eq_shock_r, i_shock_r] <- 1.0
  
  # Non-rational expectation for output gap
  # x_{t+1}^e = \lambda x_t + (1-\lambda) E_t x_{t+1}
  Gamma0[eq_xe_alt, i_var_xe_alt] <- 1.0
  Gamma0[eq_xe_alt, i_var_x] <- -1.0*lambda
  Gamma0[eq_xe_alt, i_var_xe] <- -1.0*(1.0 - lambda)
  
  # Non-rational expectation for inflation
  # \pi_{t+1}^e = \lambda \pi_t + (1-\lambda) E_t \pi_{t+1}
  Gamma0[eq_pie_alt, i_var_pie_alt] <- 1
  Gamma0[eq_pie_alt, i_var_pi] <- -1.0*lambda
  Gamma0[eq_pie_alt, i_var_pie] <- -1.0*(1.0 - lambda)
  
  # Structural shocks covariance matrix
  Sigma[i_shock_x, i_shock_x] = sigma_x * sigma_x
  Sigma[i_shock_pi, i_shock_pi] = sigma_pi * sigma_pi
  Sigma[i_shock_r, i_shock_r] = sigma_r * sigma_r
  
  # Standard Deviation of Shocks
  shocks_stddev <- sqrt(diag(Sigma))
  
  nksys.list <- list(Gamma0=Gamma0, Gamma1=Gamma1, C=C, Psi=Psi, Pi=Pi, Sigma=Sigma, shocks_stddev=shocks_stddev)
  return(nksys.list)
}

check_nksys <- function(nksys.list) {
  n <- nrow(nksys.list$Gamma0)
  nexp <- ncol(nksys.list$Pi)
  
  A <- matrix(as.complex(nksys.list$Gamma0), nrow=n, ncol=n)
  B <- matrix(as.complex(nksys.list$Gamma1), nrow=n, ncol=n)
  
  nkqz <- qz(A,B)
  eigv <- nkqz$BETA / nkqz$ALPHA
  eigv <- abs(Mod(eigv))
  
  nexplosive <- sum(abs(eigv)>=1.0)
  explosive_eigs <- abs(eigv)>=1.0
  nkqz.ord <- qz(A,B, !explosive_eigs)
  
  # Test it!
  #eigv <- nkqz.ord$BETA / nkqz.ord$ALPHA
  #eigv <- abs(Mod(eigv))
  #(eigv)
  
  retval <- "Error"
  if(nexplosive<nexp) {
    return("Indeterminacy")
  }
  
  pi_tilde <- t(nkqz.ord$Q) %*% nksys.list$Pi
  pi_tilde_2 <- pi_tilde[(n-nexplosive+1):n, ]
  rnk <- rankMatrix(pi_tilde_2)[1]
  
  if(rnk<nexp) retval <- "Indeterminacy"
  if(rnk==nexp) retval <- "Unique"
  if(rnk>nexp) retval <- "No solution"
  return(retval)
}

nkirf <- function(nksys.list, nirf=12) {
  # Solve System
  nksys <- new(gensys)
  nksys$build(Gamma0, Gamma1, C, Psi, Pi)
  nksys$solve()
  #nksys$G_sol
  #nksys$impact_sol
  
  irf.ts <- gensys_irf_ts(nksys, shocks=shocks_stdev, nirf=nirf, varnames=varnames, shocknames=shocknames)
  return(irf.ts)
}