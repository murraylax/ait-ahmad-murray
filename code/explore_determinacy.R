source("nksys.R")
library(latex2exp)

# Parameters
pistar <- (1.0 + 0.02)^{0.25} - 1.0 # 2% steady state annual inflation
rstar <- (1.0 + 0.04)^{0.25} - 1.0 # 4% steady state annual nominal interest rate

beta <- 1.0 / (1.0 + rstar) # Discount rate
sigma <- 1.0/1.38 # Preference parameter: # Schmitt-Grohe and Uribe use 1/2.0, Smets and Wouters find 1/1.38

calvo <- 0.66 # Schmitt-Grohe and Uribe use 0.8, Smets and Wouters find 0.66
kappa <- (1.0-calvo)*(1.0-calvo*beta)/calvo # Phillips curve coefficient

rho_r <- 0.0 # Interest rate smoothing parameter
psi_pi <- 1.50 # Taylor rule response to inflation target
psi_x <- 0.5 # Taylor rule response to output gap

rho_x <- 0.0 # Persistence in the demand shock
rho_pi <- 0.0 # Persistence in the supply shock

sigma_pi <- 0.25 # Shock to inflation (0.25% shock to inflation)
sigma_x <- 0.25 # Demand shock (0.25% shock to output gap)
sigma_r <- 0.25 # Shock to monetary policy (25 basis points)

# Test out the model

deltaF = 0
lambda = 0 # lambda = 0 means fully rational expectations
gamma = 0.5 # Gamma is backward-looking weight, gamma=0 means purely forward-looking
deltaB = 1.0 # Weight on most recent observation (current value)

# Test it out
#nksys.list <- nksys(gamma,deltaB,deltaF,lambda)
#solution <- check_nksys(nksys.list)
#(solution)

nsim_deltaF <- 301
deltaFsim <- seq(0.0,1.0, by=1/(nsim_deltaF-1))

nsim_gamma <- 301
gammaSim <- seq(0.0,1.0, by=1/(nsim_gamma-1))
  
all.df <- tibble(deltaF=as.double(NA), 
                 gamma=as.double(NA),
                 solution=as.character(NA),
                 unique=as.logical(NA), .rows=nsim_gamma*nsim_deltaF)

s <- 0
for(deltaF in deltaFsim) {
  ps <- sprintf("deltaF=%s\n", deltaF)
  print(ps)
  for(gamma in gammaSim) {
    s <- s + 1
    nksys.list <- nksys(gamma,deltaB,deltaF,lambda)
    solution <- check_nksys(nksys.list)
    all.df$gamma[s] <- gamma
    all.df$deltaF[s] <- deltaF
    all.df$solution[s] <- solution
    all.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

ggplot(all.df, aes(x=gamma, y=deltaF, color=solution)) +
  geom_point(alpha=0.3,stroke=0, shape=19, size=2) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=c("firebrick4", "dodgerblue4")) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  labs(title="Regions of Determinacy over Backward-Looking Weight and\nLength of Forward-Looking Windows",
       x=TeX("Backward-Window Weight: $\\gamma$"), y=TeX("Forward Weight: $\\delta_F$"),
       col="") +
  coord_cartesian(ylim=c(0,1)) +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) ->
  gg.gf
show(gg.gf)
ggsave(filename="./gamma_deltaF.png", plot=gg.gf, width=10, height=8)



# Search over lambda
# Default values
deltaF = 0
lambda = 0 # lambda = 0 means fully rational expectations
gamma = 0.0 # Gamma is backward-looking weight, gamma=0 means purely forward-looking
deltaB = 1.0 # Weight on most recent observation (current value)

nsim_deltaF <- 301
deltaFsim <- seq(0.0,1.0, by=1/(nsim_deltaF-1))

nsim_lambda <- 301
lambdaSim <- seq(0.0,1.0, by=1/(nsim_lambda-1))

all.df <- tibble(deltaF=as.double(NA), 
                 lambda=as.double(NA),
                 solution=as.character(NA),
                 unique=as.logical(NA), .rows=nsim_gamma*nsim_deltaF)

s <- 0
for(deltaF in deltaFsim) {
  ps <- sprintf("deltaF=%s\n", deltaF)
  print(ps)
  for(lambda in lambdaSim) {
    s <- s + 1
    nksys.list <- nksys(gamma,deltaB,deltaF,lambda)
    solution <- check_nksys(nksys.list)
    all.df$lambda[s] <- lambda
    all.df$deltaF[s] <- deltaF
    all.df$solution[s] <- solution
    all.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

ggplot(all.df, aes(x=lambda, y=deltaF, color=solution)) +
  geom_point(alpha=0.3,stroke=0, shape=19, size=2) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=c("firebrick4", "dodgerblue4")) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  labs(title="Regions of Determinacy for Forward Windows with Naive Expectations", 
       x=TeX("Proportion of Naive Expectations: $\\lambda$"), y=TeX("Forward Weight: $\\delta_F$"),
       x="lambda", y="delta_F", 
       col="") +
  coord_cartesian(ylim=c(0,1)) +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) ->
  gg.lf
show(gg.lf)
ggsave(filename="./lambda_deltaF.png", plot=gg.lf, width=10, height=8)


# Search over psi_pi
# Default values
deltaF = 0
lambda = 0 # lambda = 0 means fully rational expectations
gamma = 0.0 # Gamma is backward-looking weight, gamma=0 means purely forward-looking
deltaB = 1.0 # Weight on most recent observation (current value)

nsim_deltaF <- 301
deltaFsim <- seq(0.0,1.0, by=1/(nsim_deltaF-1))

nsim_psipi <- 601
minval <- 0.9
maxval <- 3.0
psipiSim <- seq(minval, maxval, by=(maxval-minval)/(nsim_psipi-1))

all.df <- tibble(deltaF=as.double(NA), 
                 psi_pi=as.double(NA),
                 solution=as.character(NA),
                 unique=as.logical(NA), .rows=nsim_psipi*nsim_deltaF)

s <- 0
for(deltaF in deltaFsim) {
  ps <- sprintf("deltaF=%s\n", deltaF)
  print(ps)
  for(psi_pi in psipiSim) {
    s <- s + 1
    nksys.list <- nksys(gamma,deltaB,deltaF,lambda,psi_pi)
    solution <- check_nksys(nksys.list)
    all.df$psi_pi[s] <- psi_pi
    all.df$deltaF[s] <- deltaF
    all.df$solution[s] <- solution
    all.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

all.df %>%
  filter(psi_pi>=0.9) %>%
  ggplot(aes(x=psi_pi, y=deltaF, color=solution)) +
  geom_point(alpha=0.3,stroke=0, shape=19, size=2) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=c("firebrick4", "dodgerblue4")) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  theme(legend.position = "bottom") +
  labs(title="Regions of Determinacy over Forward Windows and\n Monetary Policy Response to Inflation", 
       x=TeX("Taylor Rule Coefficient on Inflation: $\\psi_\\pi$"), y=TeX("Forward Weight: $\\delta_F$"),
       col="") +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) ->
  gg.pif
show(gg.pif)
ggsave(filename="./pi_deltaF.png", plot=gg.pif, width=10, height=8)


# Search over deltaB
# Default values
deltaF = 0
lambda = 0 # lambda = 0 means fully rational expectations
gamma = 0.5 # Gamma is backward-looking weight, gamma=0 means purely forward-looking
deltaB = 0.0 # Weight on most recent observation (current value)

nsim_deltaF <- 300
deltaFsim <- seq(0.001,0.999, by=1/(nsim_deltaF-1))

nsim_deltaB <- 300
deltaBsim <- seq(0.001,0.999, by=1/(nsim_deltaB-1))

all.df <- tibble(deltaF=as.double(NA), 
                 deltaB=as.double(NA),
                 solution=as.character(NA),
                 unique=as.logical(NA), .rows=nsim_deltaB*nsim_deltaF)

s <- 0
for(deltaF in deltaFsim) {
  ps <- sprintf("deltaF=%s\n", deltaF)
  print(ps)
  for(deltaB in deltaBsim) {
    s <- s + 1
    nksys.list <- nksys(gamma,deltaB,deltaF,lambda)
    solution <- check_nksys(nksys.list)
    all.df$deltaB[s] <- deltaB
    all.df$deltaF[s] <- deltaF
    all.df$solution[s] <- solution
    all.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

all.df %>%
  na.omit() %>%
  ggplot(aes(x=deltaB, y=deltaF, color=solution)) +
  geom_point(alpha=0.4,stroke=0, shape=19, size=2) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_x_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=c("firebrick4", "dodgerblue4")) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  theme(legend.position = "bottom") +
  labs(title="Regions of Determinacy for Forward & Backward Windows", x=TeX("Backward Weight: $\\delta_B$"), y=TeX("Forward Weight: $\\delta_F$"), col="") +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) +
  theme(legend.position = "bottom") ->
  gg.bf

show(gg.bf)
ggsave(filename="./deltab_deltaF.png", plot=gg.bf, width=10, height=8)


  