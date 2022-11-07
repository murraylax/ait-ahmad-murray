library(tidyverse)
library(BMR)
#devtools::install_github("kthohr/BMR")

# Put New Keynesian into Sims's GenSys form: Gamma_O y_t = C + Gamma_1 y_{t-1} + \Psi \epsilon_t + \Pi \eta_t
# where y_t = E_{t-1} y_{t-1} + \eta_t

# Indices into matrix y_t
i_x <- 1    # Output gap x_t
i_pi <- 2   # Inflation rate \pi_t
i_r <- 3    # Federal funds rate r_t
i_xe <- 4   # Expectation E_t x_{t+1}
i_pie <- 5  # Expectation E_t \pi_{t+1}
nvar <- 5

varnames <- rep(NA, nvar)
varnames[i_x] <- "Output Gap"
varnames[i_pi] <- "Inflation"
varnames[i_r] <- "Federal Funds Rate"
varnames[i_xe] <- "Expected Output Gap"
varnames[i_pie] <- "Expected Inflation Rate"

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
neq <- 5

# Parameters
pistar <- (1.0 + 0.02)^{0.25} - 1.0 # 2% steady state annual inflation
rstar <- (1.0 + 0.04)^{0.25} - 1.0 # 4% steady state annual nominal interest rate
beta <- 1.0 / (1.0 + rstar)
sigma <- 2.0
kappa <- 0.05
rho <- 0.8
psi_pi <- 1.50
psi_x <- 0.5
sigma_pi <- 0.25
sigma_x <- 0.25
sigma_r <- 0.25


# Setup matrices
Gamma0 <- matrix(nrow=neq, ncol=nvar, data=0)
C <- matrix(nrow=neq, ncol=1, data=0)
Gamma1 <- matrix(nrow=neq, ncol=nvar, data=0)
Psi <- matrix(nrow=neq, ncol=nshocks, data=0)
Pi <- matrix(nrow=neq, ncol=nexp, data=0)
Sigma <- matrix(nrow=nshocks, ncol=nshocks, data=0)


# IS equation
# x_t = E_t x_{t+1} - 1/sigma \left[(r_t - r*) - (E_t \pi_{t+1} - \pi*)\right] + \epsilon_t^x
Gamma0[eq_IS, i_x] <- 1.0
Gamma0[eq_IS, i_r] <- 1.0/sigma
Gamma0[eq_IS, i_xe] <- -1.0 
Gamma0[eq_IS, i_pie] <- -1.0/sigma
C[eq_IS, 1] <- 1.0/sigma*(rstar-pistar)
Psi[eq_IS, i_shock_x] <- 1.0

# Phillips curve
# (\pi_t - \pi*) = \beta E_t(\pi_{t+1} - \pi*) + \kappa x_t + \epsilon_t^\pi
Gamma0[eq_Phillips, i_x] <- -1.0*kappa
Gamma0[eq_Phillips, i_pi] <- 1.0
Gamma0[eq_Phillips, i_pie] <- -1.0*beta
C[eq_Phillips, 1] <- (1.0-beta)*pistar
Psi[eq_Phillips, i_shock_pi] <- 1.0

# Taylor Rule
# r_t - r* = \rho (r_{t-1} - r*) + (1-\rho) \left[ \psi_\pi (\pi_t - \pi^*) + \psi_x x_t \right] + \epsilon_t^r
Gamma0[eq_Taylor, i_x] <- -1.0*(1.0-rho)*psi_x
Gamma0[eq_Taylor, i_pi] <- -1.0*(1.0-rho)*psi_pi
Gamma0[eq_Taylor, i_r] <- 1.0
C[eq_Taylor, 1] <- (1.0-rho)*(rstar - psi_pi*pistar)
Gamma1[eq_Taylor, i_r] <- rho
Psi[eq_Taylor, i_shock_r] <- 1.0

# Output Gap Expectation
Gamma0[eq_xe, i_x] <- 1.0
Gamma1[eq_xe, i_xe] <- 1.0
Pi[eq_xe, i_exp_xe] <- 1.0


# Inflation Expectation
Gamma0[eq_pie, i_pi] <- 1.0
Gamma1[eq_pie, i_pie] <- 1.0
Pi[eq_pie, i_exp_pie] <- 1.0

# Variance of Shocks
Sigma[i_shock_x, i_shock_x] <- sigma_x*sigma_x
Sigma[i_shock_pi, i_shock_pi] <- sigma_pi*sigma_pi
Sigma[i_shock_r, i_shock_r] <- sigma_r*sigma_r

shocks_stdev <- vector(length=nshocks)
shocks_stdev[i_shock_x] <- sigma_x
shocks_stdev[i_shock_pi] <- sigma_pi
shocks_stdev[i_shock_r] <- sigma_r


nksys <- new(gensys)
nksys$build(Gamma0, Gamma1, C, Psi, Pi)
nksys$solve()
nksys$G_sol
nksys$impact_sol

nksys$shocks_cov <- Sigma

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

  