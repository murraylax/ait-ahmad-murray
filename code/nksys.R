library(tidyverse)
library(BMR)
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
nvar <- 8

varnames <- rep(NA, nvar)
varnames[i_var_x] <- "Output Gap"
varnames[i_var_pi] <- "Inflation"
varnames[i_var_r] <- "Federal Funds Rate"
varnames[i_var_xe] <- "Expected Output Gap"
varnames[i_var_pie] <- "Expected Inflation Rate"
varnames[i_var_shock_x] <- "Demand Shock"
varnames[i_var_shock_pi] <- "Cost Shock"
varnames[i_var_shock_r] <- "Monetary Policy Rate"

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
nexp <- 2

# Equation indices
eq_IS <- 1
eq_Phillips <- 2
eq_Taylor <- 3
eq_xe <- 4
eq_pie <- 5
eq_shock_x <- 6
eq_shock_pi <- 7
eq_shock_r <- 8
neq <- 8

# Parameters
pistar <- (1.0 + 0.02)^{0.25} - 1.0 # 2% steady state annual inflation
rstar <- (1.0 + 0.04)^{0.25} - 1.0 # 4% steady state annual nominal interest rate
beta <- 1.0 / (1.0 + rstar)
sigma <- 2.0
kappa <- 0.05
rho_r <- 0.8
psi_pi <- 1.50
psi_x <- 0.5
sigma_pi <- 0.25
sigma_x <- 0.25
sigma_r <- 0.25
rho_x <- 0.7
rho_pi <- 0.7

# Setup matrices
Gamma0 <- matrix(nrow=neq, ncol=nvar, data=0)
C <- matrix(nrow=neq, ncol=1, data=0)
Gamma1 <- matrix(nrow=neq, ncol=nvar, data=0)
Psi <- matrix(nrow=neq, ncol=nshocks, data=0)
Pi <- matrix(nrow=neq, ncol=nexp, data=0)
Sigma <- matrix(nrow=nshocks, ncol=nshocks, data=0)


# IS equation
# x_t = E_t x_{t+1} - 1/sigma \left[(r_t - r*) - (E_t \pi_{t+1} - \pi*)\right] + \epsilon_t^x
Gamma0[eq_IS, i_var_x] <- 1.0
Gamma0[eq_IS, i_var_r] <- 1.0/sigma
Gamma0[eq_IS, i_var_xe] <- -1.0 
Gamma0[eq_IS, i_var_pie] <- -1.0/sigma
C[eq_IS, 1] <- 1.0/sigma*(rstar-pistar)
Gamma0[eq_IS, i_var_shock_x] <- -1.0

# Phillips curve
# (\pi_t - \pi*) = \beta E_t(\pi_{t+1} - \pi*) + \kappa x_t + \epsilon_t^\pi
Gamma0[eq_Phillips, i_var_x] <- -1.0*kappa
Gamma0[eq_Phillips, i_var_pi] <- 1.0
Gamma0[eq_Phillips, i_var_pie] <- -1.0*beta
C[eq_Phillips, 1] <- (1.0-beta)*pistar
Gamma0[eq_Phillips, i_var_shock_pi] <- -1.0

# Taylor Rule
# r_t - r* = \rho (r_{t-1} - r*) + (1-\rho) \left[ \psi_\pi (\pi_t - \pi^*) + \psi_x x_t \right] + \epsilon_t^r
Gamma0[eq_Taylor, i_var_x] <- -1.0*(1.0-rho)*psi_x
Gamma0[eq_Taylor, i_var_pi] <- -1.0*(1.0-rho)*psi_pi
Gamma0[eq_Taylor, i_var_r] <- 1.0
C[eq_Taylor, 1] <- (1.0-rho)*(rstar - psi_pi*pistar)
Gamma1[eq_Taylor, i_var_r] <- rho
Gamma0[eq_Taylor, i_var_shock_r] <- -1.0

# Output Gap Expectation
Gamma0[eq_xe, i_var_x] <- 1.0
Gamma1[eq_xe, i_var_xe] <- 1.0
Pi[eq_xe, i_exp_xe] <- 1.0

# Inflation Expectation
Gamma0[eq_pie, i_var_pi] <- 1.0
Gamma1[eq_pie, i_var_pie] <- 1.0
Pi[eq_pie, i_exp_pie] <- 1.0

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

# Standard Deviation of Shocks
shocks_stdev <- vector(length=nshocks)
shocks_stdev[i_shock_x] <- sigma_x
shocks_stdev[i_shock_pi] <- sigma_pi
shocks_stdev[i_shock_r] <- sigma_r

# Solve System
nksys <- new(gensys)
nksys$build(Gamma0, Gamma1, C, Psi, Pi)
nksys$solve()
nksys$G_sol
nksys$impact_sol

nirf <- 12
irf.ts <- gensys_irf_ts(nksys, shocks=shocks_stdev, nirf=nirf, varnames=varnames, shocknames=shocknames)

irf.ts %>%
  filter(Variable %in% c("Output Gap", "Inflation", "Federal Funds Rate")) %>%
  ggplot(aes(x=Period, y=Impact)) +
  geom_line(size=1.2, color="dodgerblue4") +
  facet_grid(Shock~Variable) +
  theme_bw() + 
  theme(text = element_text(size=16)) +
  theme(strip.background = element_rect(fill="white")) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(label=percent_format(scale=1, accuracy=0.1))

  