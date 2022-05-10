source("nksys.R")
library(latex2exp)

# Plot parameters
shape <- 19 # Circle
size <- 2.5 # Size of the point
#scale_colors <- c("firebrick4", "dodgerblue4")
scale_colors <- c("firebrick2", "dodgerblue2")
scale_colors <- c("gray47", "gray87")
alpha <- 1


# Parameters
pistar <- (1.0 + 0.02)^{0.25} - 1.0 # 2% steady state annual inflation
rstar <- (1.0 + 0.04)^{0.25} - 1.0 # 4% steady state annual nominal interest rate

beta <- 1.0 / (1.0 + rstar) # Discount rate
sigma <- 1.0/1.38 # Preference parameter: # Schmitt-Grohe and Uribe use 1/2.0, Smets and Wouters find 1/1.38

calvo <- 0.66 # Schmitt-Grohe and Uribe use 0.8, Smets and Wouters find 0.66
kappa <- (1.0-calvo)*(1.0-calvo*beta)/calvo # Phillips curve coefficient

#rho_r <- 0.81 # Interest rate smoothing parameter, Smets and Wouters
rho_r <- 0.0 # Interest rate smoothing parameter, Smets and Wouters
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

nsim_deltaF <- 201
minval <- 0
maxval <- 0.5
deltaFsim <- seq(minval, maxval, by=(maxval-minval)/(nsim_deltaF-1))

nsim_gamma <- 201
minval <- 0
maxval <- 0.99
gammaSim <- seq(minval, maxval, by=(maxval-minval)/(nsim_gamma-1))
  
gamma_deltaF.df <- tibble(deltaF=as.double(NA), 
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
    solution <- checksys(nksys.list)
    gamma_deltaF.df$gamma[s] <- gamma
    gamma_deltaF.df$deltaF[s] <- deltaF
    gamma_deltaF.df$solution[s] <- solution
    gamma_deltaF.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

ggplot(gamma_deltaF.df, aes(x=gamma, y=deltaF, color=solution)) +
  geom_point(alpha=alpha, stroke=0, shape=shape, size=size) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=scale_colors) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  theme(legend.position="bottom") +
  labs(title="Regions of Determinacy over Backward-Looking Weight and\nLength of Forward-Looking Windows",
       x=TeX("Backward-Window Weight: $\\gamma$"), y=TeX("Forward Weight: $\\delta_F$"),
       col="") +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) ->
  gg.gf
show(gg.gf)
ggsave(filename="./gamma_deltaF.png", plot=gg.gf, width=10, height=8)

gg.gf.notitle <- gg.gf + labs(title="") + theme(legend.position = "none")
show(gg.gf.notitle)
ggsave(filename="./gamma_deltaF_notitle.png", plot=gg.gf.notitle)



# Search over lambda
# Default values
deltaF = 0
lambda = 0 # lambda = 0 means fully rational expectations
gamma = 0.0 # Gamma is backward-looking weight, gamma=0 means purely forward-looking
deltaB = 1.0 # Weight on most recent observation (current value)

nsim_deltaF <- 201
delta_low <- 0.01
delta_high <- 0.5
deltaFsim <- seq(delta_low, delta_high, by=(delta_high-delta_low)/(nsim_deltaF-1))

nsim_lambda <- 201
lambda_low <- 0.01
lambda_high <- 0.99
lambdaSim <- seq(lambda_low,lambda_high, by=(lambda_high-lambda_low)/(nsim_lambda-1))

lambda_deltaF.df <- tibble(deltaF=as.double(NA), 
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
    solution <- checksys(nksys.list)
    lambda_deltaF.df$lambda[s] <- lambda
    lambda_deltaF.df$deltaF[s] <- deltaF
    lambda_deltaF.df$solution[s] <- solution
    lambda_deltaF.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

ggplot(lambda_deltaF.df, aes(x=lambda, y=deltaF, color=solution)) +
  geom_point(alpha=alpha, stroke=0, shape=shape, size=size) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=scale_colors) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  theme(legend.position="bottom") +
  labs(title="Regions of Determinacy for Forward Windows with Naive Expectations", 
       x=TeX("Proportion of Naive Expectations: $\\lambda$"), y=TeX("Forward Weight: $\\delta_F$"),
       x="lambda", y="delta_F", 
       col="") +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) ->
  gg.lf
show(gg.lf)
ggsave(filename="./lambda_deltaF.png", plot=gg.lf, width=10, height=8)

gg.lf.notitle <- gg.lf + labs(title="") + theme(legend.position = "none")
show(gg.lf.notitle)
ggsave(filename="./lambda_deltaF_notitle.png", plot=gg.lf.notitle)



# Search over psi_pi
# Default values
deltaF = 0
lambda = 0 # lambda = 0 means fully rational expectations
gamma = 0.0 # Gamma is backward-looking weight, gamma=0 means purely forward-looking
deltaB = 1.0 # Weight on most recent observation (current value)

nsim_deltaF <- 201
minval <- 0.01
maxval <- 0.5
deltaFsim <- seq(minval,maxval, by=(maxval-minval)/(nsim_deltaF-1))

nsim_psipi <- 201
minval <- 0.9
maxval <- 3.0
psipiSim <- seq(minval, maxval, by=(maxval-minval)/(nsim_psipi-1))

psipi_deltaF.df <- tibble(deltaF=as.double(NA), 
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
    solution <- checksys(nksys.list)
    psipi_deltaF.df$psi_pi[s] <- psi_pi
    psipi_deltaF.df$deltaF[s] <- deltaF
    psipi_deltaF.df$solution[s] <- solution
    psipi_deltaF.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

ggplot(psipi_deltaF.df, aes(x=psi_pi, y=deltaF, color=solution)) +
  geom_point(alpha=alpha, stroke=0, shape=shape, size=size) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=scale_colors) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  theme(legend.position = "bottom") +
  labs(title="Regions of Determinacy over Forward Windows and\n Monetary Policy Response to Inflation", 
       x=TeX("Taylor Rule Coefficient on Average Inflation: $\\psi_\\pi$"), y=TeX("Forward Weight: $\\delta_F$"),
       col="") +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) ->
  gg.pif
show(gg.pif)
ggsave(filename="./pi_deltaF.png", plot=gg.pif, width=10, height=8)

gg.pif.notitle <- gg.pif + labs(title="") + theme(legend.position = "none")
show(gg.pif.notitle)
ggsave(filename="./pif_deltaF_notitle.png", plot=gg.pif.notitle)


# Search over deltaB
# Default values
deltaF = 0
lambda = 0 # lambda = 0 means fully rational expectations
gamma = 0.25 # Gamma is backward-looking weight, gamma=0 means purely forward-looking
deltaB = 0.0 # Weight on most recent observation (current value)

nsim_deltaF <- 201
minval <- 0.01
maxval <- 0.5
deltaFsim <- seq(minval, maxval, by=(maxval-minval)/(nsim_deltaF-1))

nsim_deltaB <- 201
minval <- 0.01
maxval <- 0.99
deltaBsim <- seq(minval, maxval, by=(maxval-minval)/(nsim_deltaB-1))

deltaB_deltaF.df <- tibble(deltaF=as.double(NA), 
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
    solution <- checksys(nksys.list)
    deltaB_deltaF.df$deltaB[s] <- deltaB
    deltaB_deltaF.df$deltaF[s] <- deltaF
    deltaB_deltaF.df$solution[s] <- solution
    deltaB_deltaF.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

ggplot(deltaB_deltaF.df, aes(x=deltaB, y=deltaF, color=solution)) +
  geom_point(alpha=alpha, stroke=0, shape=shape, size=size) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_x_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=scale_colors) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  theme(legend.position = "bottom") +
  labs(title="Regions of Determinacy for Forward & Backward Windows", x=TeX("Backward Weight: $\\delta_B$"), y=TeX("Forward Weight: $\\delta_F$"), col="") +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) +
  theme(legend.position = "bottom") ->
  gg.bf

show(gg.bf)
ggsave(filename="./deltab_deltaF.png", plot=gg.bf, width=10, height=8)

gg.bf.notitle <- gg.bf + labs(title="") + theme(legend.position = "none")
show(gg.bf.notitle)
ggsave(filename="./deltab_deltaF_notitle.png", plot=gg.bf.notitle)



# Search over rho_r
# Default values
deltaF = 0
lambda = 0 # lambda = 0 means fully rational expectations
gamma = 0.0 # Gamma is backward-looking weight, gamma=0 means purely forward-looking
deltaB = 0.0 # Weight on most recent observation (current value)

nsim_rho <- 201
minval <- 0.01
maxval <- 0.99
rhosim <- seq(minval, maxval, by=(maxval-minval)/(nsim_rho-1))

nsim_deltaF <- 201
minval <- 0.01
maxval <- 0.5
deltaFsim <- seq(minval, maxval, by=(maxval-minval)/(nsim_deltaF-1))

rho_deltaF.df <- tibble(deltaF=as.double(NA), 
                           rho=as.double(NA),
                           solution=as.character(NA),
                           unique=as.logical(NA), .rows=nsim_rho*nsim_deltaF)

s <- 0
for(deltaF in deltaFsim) {
  ps <- sprintf("deltaF=%s\n", deltaF)
  print(ps)
  for(rho in rhosim) {
    s <- s + 1
    nksys.list <- nksys(gamma,deltaB,deltaF,lambda,rho_r=rho)
    solution <- checksys(nksys.list)
    rho_deltaF.df$rho[s] <- rho
    rho_deltaF.df$deltaF[s] <- deltaF
    rho_deltaF.df$solution[s] <- solution
    rho_deltaF.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

ggplot(rho_deltaF.df, aes(x=rho, y=deltaF, color=solution)) +
  geom_point(alpha=alpha, stroke=0, shape=shape, size=size) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_x_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=scale_colors) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  theme(legend.position = "bottom") +
  labs(title="Regions of Determinacy for Forward Windows and \nMonetary Policy Persistence", x=TeX("Monetary Policy Persistence: $\\rho_r$"), y=TeX("Forward Weight: $\\delta_F$"), col="") +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) +
  theme(legend.position = "bottom") ->
  gg.rf

show(gg.rf)
ggsave(filename="./rho_deltaF.png", plot=gg.rf, width=10, height=8)

gg.rf.notitle <- gg.rf + labs(title="") + theme(legend.position = "none")
show(gg.rf.notitle)
ggsave(filename="./rho_deltaF_notitle.png", plot=gg.rf.notitle)






# Search over psi_x
# Default values
deltaF = 0
lambda = 0 # lambda = 0 means fully rational expectations
gamma = 0.0 # Gamma is backward-looking weight, gamma=0 means purely forward-looking
deltaB = 0.0 # Weight on most recent observation (current value)

nsim_psix <- 201
minval <- 0.01
maxval <- 0.99
psi_x_sim <- seq(minval, maxval, by=(maxval-minval)/(nsim_psix-1))

nsim_deltaF <- 201
minval <- 0.01
maxval <- 0.5
deltaFsim <- seq(minval, maxval, by=(maxval-minval)/(nsim_deltaF-1))

psix_deltaF.df <- tibble(deltaF=as.double(NA), 
                        psi_x=as.double(NA),
                        solution=as.character(NA),
                        unique=as.logical(NA), .rows=nsim_rho*nsim_deltaF)

s <- 0
for(deltaF in deltaFsim) {
  ps <- sprintf("deltaF=%s\n", deltaF)
  print(ps)
  for(psi_x in psi_x_sim) {
    s <- s + 1
    nksys.list <- nksys(gamma,deltaB,deltaF,lambda,psi_x=psi_x)
    solution <- checksys(nksys.list)
    psix_deltaF.df$psi_x[s] <- psi_x
    psix_deltaF.df$deltaF[s] <- deltaF
    psix_deltaF.df$solution[s] <- solution
    psix_deltaF.df$unique[s] <- solution=="Unique"
  }
}

# Plot determinacy / indeterminacy region

ggplot(psix_deltaF.df, aes(x=psi_x, y=deltaF, color=solution)) +
  geom_point(alpha=alpha, stroke=0, shape=shape, size=size) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_x_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=scale_colors) +
  theme_bw() +
  theme(text=element_text(size=18)) +
  theme(legend.position = "bottom") +
  labs(title="Regions of Determinacy for Forward Windows and \nMonetary Policy Response to Output Gap", x=TeX("Monetary Policy Response to Output: $\\psi_x$"), y=TeX("Forward Weight: $\\delta_F$"), col="") +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) +
  theme(legend.position = "bottom") ->
  gg.psixf

show(gg.psixf)
ggsave(filename="./psix_deltaF.png", plot=gg.psixf, width=10, height=8)

gg.psixf.notitle <- gg.psixf + labs(title="") + theme(legend.position = "none")
show(gg.psixf.notitle)
ggsave(filename="./psix_deltaF_notitle.png", plot=gg.psixf.notitle)




## Arrange grid ------------------------------------------------
library(gridExtra)
grid.arrange(gg.bf.notitle, gg.gf.notitle, gg.lf.notitle, gg.pif.notitle)

grid.arrange(gg.bf, gg.gf, gg.lf, gg.pif)
## -------------------------------------------------------------

df <- bind_rows(deltaB_deltaF.df, gamma_deltaF.df)
df <- bind_rows(df, lambda_deltaF.df)
df <- bind_rows(df, psipi_deltaF.df)
df <- bind_rows(df, psix_deltaF.df)
df <- bind_rows(df, rho_deltaF.df)

df.gather <- gather(df, key="param", value="value", c(deltaB, gamma, lambda, psi_pi, psi_x, rho))

parmlabels <- c(TeX("Panel (B): Backward Window, $\\delta_B$"), 
                TeX("Panel (A): Backward Weight, $\\gamma$"), 
                TeX("Panel (C): Proportion Naive, $\\lambda$"), 
                TeX("Panel (D): Taylor Rule Inflation, $\\psi_\\pi$"),
                TeX("Panel (E): Taylor Rule Output Gap, $\\psi_x$"),
                TeX("Panel (F): Taylor Rule Persistence, $\\rho_r$")
)

df.gather$Parameter <- as.factor(df.gather$param)
levels(df.gather$Parameter) <- parmlabels

sort_parmlabels <- c(
                TeX("Panel (A): Backward Weight, $\\gamma$"), 
                TeX("Panel (B): Backward Window, $\\delta_B$"), 
                TeX("Panel (C): Proportion Naive, $\\lambda$"), 
                TeX("Panel (D): Taylor Rule Inflation, $\\psi_\\pi$"),
                TeX("Panel (E): Taylor Rule Output Gap, $\\psi_x$"),
                TeX("Panel (F): Taylor Rule Persistence, $\\rho_r$")
                )

df.gather$Parameter <- factor(df.gather$Parameter, levels=sort_parmlabels, ordered=TRUE)



ggplot(df.gather, aes(x=value, y=deltaF, color=solution)) +
  geom_point(alpha=alpha, stroke=0, shape=shape, size=size) +
  scale_y_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_x_continuous(label=number_format(accuracy=0.01), breaks=pretty_breaks(n=5)) +
  scale_color_manual(values=scale_colors) +
  facet_wrap(~Parameter, nrow=2, scales="free_x", labeller = "label_parsed", strip.position="bottom") +
  theme_bw() +
  theme(text=element_text(size=14)) +
  theme(legend.position = c(0.5,-0.18), legend.direction = "horizontal") +
  theme(plot.margin = margin(b = 22, t=5, r=5, l=5)) +
  theme(panel.spacing = unit(1, "lines")) +
  theme(strip.placement = "outside", strip.background = element_rect(size=0, fill="white")) +
  labs(title="Regions of Determinacy for Forward-Looking Windows", x="", y=TeX("Forward Weight: $\\delta_F$"), col="") +
  guides(color = guide_legend(override.aes = list(size=10, alpha=1))) ->
  gg.all

show(gg.all)
ggsave(filename="./determinacy.png", plot=gg.all, width=10, height=6.5)

gg.all.notitle <- gg.all + labs(title="") + theme(plot.margin=margin(b = 25, t=-10, r=5, l=5))
show(gg.all.notitle)
ggsave(filename="./determinacy_notitle.png", plot=gg.all.notitle, width=10, height=6.5)

# Compute maximum values of indeterminacy
df.gather %>%
  group_by(param) %>%
  filter(!is.na(value)) %>%
  filter(!(param=="psi_pi" & value<1.0)) %>%
  filter(solution=="Indeterminacy") %>%
  filter(deltaF==max(deltaF)) %>%
  select(param, deltaF, value) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(value=round(value,1)) %>%
  mutate(nquarters = 1/deltaF) %>% 
  select(param, value, deltaF, nquarters) ->
  tb

tb


