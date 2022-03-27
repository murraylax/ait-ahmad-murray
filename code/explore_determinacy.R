source("nksys.R")

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

# For every gamma, find max deltaF where unique==TRUE
all.df %>%
  group_by(gamma) %>%
  filter(unique) %>%
  summarise(minDeltaF = min(deltaF)) ->
  plot.df

plot.df.sub <- filter(plot.df, gamma>0, minDeltaF>0)
plot.lm <- lm(minDeltaF ~ gamma, plot.df.sub)
plot.df$minDeltaF_smooth <- plot.df$minDeltaF
plot.df$minDeltaF_smooth[ plot.df$minDeltaF>0 & plot.df$gamma>0 ] <- plot.lm$fitted.values

ggplot(plot.df, aes(x=gamma)) +
  geom_line(mapping=aes(y=minDeltaF_smooth), color="dodgerblue4", size=1) +
  geom_ribbon(mapping=aes(x=gamma, ymin=minDeltaF_smooth, ymax=1.0), fill="dodgerblue4", alpha=0.3) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  labs(title="Regions of Determinacy for Forward Windows", x="gamma", y="delta_F") +
  coord_cartesian(ylim=c(0,1))
