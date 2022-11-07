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
i_var_piA <- 9      # Inflation target, i.e. some average
i_var_piB <- 10     # Backward-looking target
i_var_piF <- 11     # Forward-looking target
i_var_piFe <- 12    # Expected Forward-looking target forecast error
nvar <- 12

varnames <- rep(NA, nvar)
varnames[i_var_x] <- "Output Gap"
varnames[i_var_pi] <- "Inflation"
varnames[i_var_r] <- "Federal Funds Rate"
varnames[i_var_xe] <- "Expected Output Gap"
varnames[i_var_pie] <- "Expected Inflation Rate"
varnames[i_var_shock_x] <- "Demand Shock"
varnames[i_var_shock_pi] <- "Cost Shock"
varnames[i_var_shock_r] <- "Monetary Policy Rate"
varnames[i_var_piA] <- "Inflation Target"
varnames[i_var_piB] <- "Backward-Looking Average Inflation"
varnames[i_var_piF] <- "Forward-Looking Average Inflation"
varnames[i_var_piFe] <- "Expected Forward-Looking Average Inflation"


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
neq <- 12

# Parameters
pistar <- (1.0 + 0.02)^{0.25} - 1.0 # 2% steady state annual inflation
rstar <- (1.0 + 0.04)^{0.25} - 1.0 # 4% steady state annual nominal interest rate
beta <- 1.0 / (1.0 + rstar) # Discount rate
sigma <- 2.0 # Preference parameter
kappa <- 0.01 # Phillips curve coefficient
rho_r <- 0.7 # Interest rate smoothing parameter
psi_pi <- 1.50 # Taylor rule response to inflation target
psi_x <- 0.5 # Taylor rule response to output gap
sigma_pi <- 0.25 # Shock to inflation (0.25% shock to inflation)
sigma_x <- 0.25 # Demand shock (0.25% shock to output gap)
sigma_r <- 0.25 # Shock to monetary policy (25 basis points)
rho_x <- 0.7 # Persistence in the demand shock
rho_pi <- 0.7 # Persistence in the supply shock

nkmodel <- function(gamma, deltaB, deltaF)
{
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
  
  return(irf.ts)
}

gamma <- 0.3 # Backward looking weight
deltaB <- 0.25 # 0.25 weight to most recent quarter, backward looking window ~ 4 quarters
deltaF <- 0.125 # 0.125 weight to next period's expectation, forward looking window ~ 8 quarters

if(FALSE) {
  irf.ts <- nkmodel(gamma, deltaB, deltaF)
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
}

# Forward looking windows
gamma <- 0.5 # Backward looking weight
deltaB <- 0.25 # 0.25 weight to most recent quarter, backward looking window ~ 4 quarters
deltaF <- 0.25 # 0.125 weight to next period's expectation, forward looking window ~ 8 quarters

irf_all.df <- tibble()

for(deltaF in c(0.0625, 0.25, 0.5, 1)) {
  irf <- nkmodel(gamma, deltaB, deltaF)
  irf <- as_tibble(irf)
  irf$`Forward Window` <- as.integer(1/deltaF)
  irf_all.df <- bind_rows(irf_all.df, irf)
}

irf_all.df$`Forward Window` <- as.factor(irf_all.df$`Forward Window`)
levs <- sprintf("%s Quarters", levels(irf_all.df$`Forward Window`))
levels(irf_all.df$`Forward Window`) <- levs
irf_all.df$`Forward Window` <- factor(irf_all.df$`Forward Window`, levels=levs, ordered=TRUE)


mycols <- viridis::plasma(n=5)

irf_all.df %>%
  filter(Variable %in% c("Output Gap", "Inflation", "Federal Funds Rate")) %>%
  ggplot(aes(x=Period, y=Impact, col=`Forward Window`)) +
  geom_line(size=1.2) +
  facet_grid(Shock~Variable) +
  theme_bw() + 
  theme(text = element_text(size=16)) +
  theme(strip.background = element_rect(fill="white")) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(label=percent_format(scale=1, accuracy=0.1)) +
  scale_color_manual(values=mycols) +
  labs(x="Periods After Shock", y="", title="Impulse Responses With Varying Forward-Looking Windows")


# Backward looking windows
gamma <- 0.5 # Backward looking weight
deltaB <- 0.25 # 0.25 weight to most recent quarter, backward looking window ~ 4 quarters
deltaF <- 0.25 # 0.25 weight to next period's expectation, forward looking window ~ 8 quarters

irf_all.df <- tibble()

for(deltaB in c(0.0625, 0.25, 0.5, 1)) {
  irf <- nkmodel(gamma, deltaB, deltaF)
  irf <- as_tibble(irf)
  irf$`Backward Window` <- as.integer(1/deltaB)
  irf_all.df <- bind_rows(irf_all.df, irf)
}

irf_all.df$`Backward Window` <- as.factor(irf_all.df$`Backward Window`)
levs <- sprintf("%s Quarters", levels(irf_all.df$`Backward Window`))
levels(irf_all.df$`Backward Window`) <- levs
irf_all.df$`Backward Window` <- factor(irf_all.df$`Backward Window`, levels=levs, ordered=TRUE)


mycols <- viridis::plasma(n=5)

irf_all.df %>%
  filter(Variable %in% c("Output Gap", "Inflation", "Federal Funds Rate")) %>%
  ggplot(aes(x=Period, y=Impact, col=`Backward Window`)) +
  geom_line(size=1.2) +
  facet_grid(Shock~Variable) +
  theme_bw() + 
  theme(text = element_text(size=16)) +
  theme(strip.background = element_rect(fill="white")) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(label=percent_format(scale=1, accuracy=0.1)) +
  scale_color_manual(values=mycols) +
  labs(x="Periods After Shock", y="", title="Impulse Responses With Varying Backward-Looking Windows")

