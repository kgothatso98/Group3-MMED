# Load required library
library(deSolve)
library(dplyr)

# Define the SIR model with mortality
sir_with_mortality <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S_w + I_w + R_w + S_o + I_o + R_o
    
    # Differential equations for whites
    dS_wdt <- -beta_ww * S_w * I_w / N_w - beta_wo * S_w * I_o/N_o
    dI_wdt <- beta_ww * S_w * I_w / N_w + beta_wo * S_w * I_o/N_o - gamma_w * I_w  - mu_w * I_w
    dR_wdt <- gamma_w * I_w
    dD_w <- mu_w * I_w  # Cumulative deaths
    
    # Differential equations for others
    dS_odt <- -beta_oo * S_o * I_o / N_o - beta_ow * S_o * I_w/N_w
    dI_odt <-  beta_oo * S_o * I_o / N_o + beta_ow * S_o * I_w/N_w- gamma_o * I_o  - mu_o * I_o
    dR_odt <- gamma_o * I_o
    dD_o <- mu_o * I_o # Cumulative deaths
    list(c( dS_wdt, dI_wdt, dR_wdt, dD_w,dS_wdt, dI_wdt, dR_wdt, dD_w))
  })
}

# Initial values
init_state <- c(S_w = 63295, I_w = 1, R_w = 0, D_w = 0, S_o = 167303, I_o = 1, R_o = 0, D_o = 0)

# Parameters
params <- c(
  beta_oo =0.6,
  beta_ow =0.2,
  beta_ww =0.6,
  beta_wo =0.2,
    # Infection rate
  gamma_o = 1/8,
  gamma_w = 1/8,
  # Recovery rate
  mu_w = 0.01 ,
  mu_o = 0.01,# Mortality rate
  
  N_o =167304,
  N_w = 63296
)

################################################################################
# Times (same as your observed data time range)
times <- seq(0, 160, by = 1)

# Your observed data (must include time and daily_deaths columns)
# obsDat <- read.csv("your_data.csv")
# Solve ODE
output <- as.data.frame(ode(y = init_state, times = times, func = sir_with_mortality, parms = params))

# Add daily_death_w, daily_deaths_o and daily deaths column (difference between successive cumulative deaths)
output <- output %>%
  mutate(
    daily_deaths_w = c(NA, diff(D_w)),
    daily_deaths_o = c(NA, diff(D_o)),
    daily_deaths_total = daily_deaths_w + daily_deaths_o
  )

# Print a preview
head(output, 10)
View(output)

###############################################################################
## Function that makes a list of disease parameters with default values
disease_params <- function(beta = 0.6 ## transmission coefficient when prevalence is 0
                           ,mu = 1/60 
                           ## 60 year natural life expectancy
){
  return(as.list(environment()))
}
#################################################################################
#Creating a function for the SIR model
SIRmod <- function(tt, yy, parms) with(c(parms,as.list(yy)), {
  ## State variables are: S, I, and R
  N <- S + I + R  
  gamma = 1/8 ## total population
  ## state variable derivatives (ODE system)
  deriv <- rep(NA,4)
  deriv[1] <- -beta*S*I/N ## Instantaneous rate of change: Susceptibles
  deriv[2] <-	beta*S*I/N - mu*I - gamma * I ## Instantaneous rate of change: Infection class I1
  deriv[3] <- gamma*I  ##Cumulative recovered people
  deriv[4] <- mu * I  # Cumulative deaths
  return(list(deriv))
})


#################################################################################
simEpidemic <- function(init_state, parms) {
  as.data.frame(ode(
    y = init_state,
    times = 0:160,   # or whatever you need
    func = sir_with_mortality,  # your new SIR function
    parms = parms
  ))
}

###########################################
##############################################
## Return the negative log of likelihood by convention


nllikelihood <- function(parms = disease_params(), obsDat = FLUepicurve[-c(1:28), ]) {
  simDat <- simEpidemic(init_state, parms = parms)
  
  # Build a data frame from sim output
  simDf <- data.frame(time = simDat$time, deaths = simDat$daily_deaths_w)
  
  # Merge by time (only keep overlapping times)
  merged <- merge(obsDat[, c("time", "ma_7day")], simDf, by = "time")
  
  # Remove rows where either observed or simulated deaths are NA
  merged <- merged[!is.na(merged$ma_7day) & !is.na(merged$deaths), ]
  
  # Optional: set negative death values to zero (if needed)
  merged$deaths[merged$deaths < 0] <- 0
  
  print(merged$deaths)
  print(merged$ma_7day)
  
  # Compute negative log-likelihood
  nlls <- suppressWarnings(-dpois(merged$ma_7day, lambda = merged$deaths, log = TRUE))
  
  
  return(sum(nlls))
}

disease_params <- function(beta = 0.5, mu = 0.01) {
  c(
    beta_oo = beta,
    beta_ow = beta / 3,
    beta_ww = beta,
    beta_wo = beta / 3,
    gamma_o = 1/8,
    gamma_w = 1/8,
    mu_o = mu,
    mu_w = mu,
    N_o = 167304,
    N_w = 63296
  )
}


nllikelihood(disease_params(beta = 3, mu = 0.01))  ## vs some random guessed parameters


#################################################################################
subsParms <- function(fit.params, fixed.params=disease_params()) {
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })
}
## Make likelihood a function of fixed and fitted parameters.
objFXN <- function(fit.params ## paramters to fit
                   , fixed.params =disease_params() ## fixed paramters
                   , obsDat=FLUepicurve[-c(1:28),]) {
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat) ## then call likelihood
}

guess.params <- c(log_beta = log(0.6), log_mu = log(1/60))
guess.params
subsParms(guess.params, disease_params())

objFXN(guess.params, disease_params())

## Select initial values for fitted parameters from which optimization routine
## will start. If you select bad initial values the algorithm can get stuck on a
## bad set of parameters. You can always try the true values as a starting point
## for this problem, although that's rarely possible in real problems.

init.pars <- c(log_beta = log(0.6), log_mu = log(1/60))
## We will start with SANN optimization since it is stochastic and therefore
## less likely to get stuck in a local minima. But then finish with Nelder-Mead
## optimization which is much faster.

###  NOTE: for trace >0 you see more progress report, bigger numbers show more
###  update
trace <- 3
trace
## SANN: This is stochastic, be CAREFUL -- sometimes it gets stuck at local minima
## for unreasonable parameters. If you see this happen, run it again!
optim.vals <- optim(par = init.pars
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = FLUepicurve[-c(1:28),]
                    , control = list(trace = trace, maxit = 150)
                    , method = "SANN")
exp(unname(optim.vals$par))

## We feed the last parameters of SANN in as the first values of Nelder-Mead
optim.vals <- optim(par = optim.vals$par
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = FLUepicurve[-c(1:28),]
                    , control = list(trace = trace, maxit = 800, reltol = 10^-7)
                    , method = "Nelder-Mead"
                    , hessian = T)
optim.vals
MLEfits <- optim.vals$par
MLEfits

exp(unname(MLEfits))

log_alpha.fit <- MLEfits["log_alpha"]
log_alpha.fit
log_Beta.fit <- MLEfits["log_Beta"]
log_Beta.fit
## Look at the output of optim. Understand what it means. Did the algorithm
## converge? Look at ?optim to understand it.

## Plot MLE fit time series
par(bty='n', lwd = 2, las = 1)
fitDat <- simEpidemic(init_state, parms = subsParms(optim.vals$par,  disease_params()))
fitDat

ggplot(fitDat, aes(x= time, y = daily_deaths)) + 
  geom_line() + 
  geom_line(data = FLUepicurve, aes(x = time, y = ma_7day), col= "red")
with(fitDat, plot(time, daily_deaths, col='blue',lty = ))
points(myDat$time, myDat$sampPrev, col = 'red', pch = 16, cex = 2)
arrows(myDat$time, myDat$uci, myDat$time, myDat$lci, col = 'red', len = .025, angle = 90, code = 3)
legend("topleft", c('truth', 'observed', 'fitted'), lty = c(1, NA, 1), pch = c(NA,16, NA),
       col = c('red', 'red', 'blue'), bty = 'n')
##################################################################################
# Initial state for the White population
init_state <- c(S = 63296, I = 1, R = 0, D = 0)
# Your observed data (must include time and daily_deaths columns)
# obsDat <- read.csv("your_data.csv")
# Solve ODE
output <- as.data.frame(ode(y = init_state, times = times, func = SIRmod, parms = params))
# Add daily deaths column (difference between successive cumulative deaths)

#################################################################################
#################################################################################
simEpidemic <- function(init = init_state, tseq = times, modFunction=SIRmod, parms = disease_params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  simDat$daily_deaths = ceiling(c(0, diff(simDat$D)))
  return(simDat)
}

## Return the negative log of likelihood by convention
nllikelihood <- function(parms = disease_params(), obsDat=FLUepicurve_W[-c(1:28),]) {
  simDat <- simEpidemic(init_state, parms=parms)
  matchedTimes <- simDat$time %in% obsDat$time
  nlls <- -dpois(obsDat$ma_7day, lambda = simDat$daily_deaths[matchedTimes], log = T)
  return(sum(nlls))
}
nllikelihood(disease_params(beta = 3, mu = 0.01))  ## vs some random guessed parameters
#################################################################################
subsParms <- function(fit.params, fixed.params=disease_params()) {
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })
}
## Make likelihood a function of fixed and fitted parameters.
objFXN <- function(fit.params ## paramters to fit
                   , fixed.params =disease_params() ## fixed paramters
                   , obsDat=FLUepiCurve_W[-c(1:28),]) {
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat) ## then call likelihood
}
objFXN(guess.params, disease_params())

## Select initial values for fitted parameters from which optimization routine
## will start. If you select bad initial values the algorithm can get stuck on a
## bad set of parameters. You can always try the true values as a starting point
## for this problem, although that's rarely possible in real problems.

init.pars <- c(log_beta = log(0.6), log_mu = log(1/60))
## We will start with SANN optimization since it is stochastic and therefore
## less likely to get stuck in a local minima. But then finish with Nelder-Mead
## optimization which is much faster.

###  NOTE: for trace >0 you see more progress report, bigger numbers show more
###  update
trace <- 3
trace
## SANN: This is stochastic, be CAREFUL -- sometimes it gets stuck at local minima
## for unreasonable parameters. If you see this happen, run it again!
optim.vals <- optim(par = init.pars
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = FLUepiCurve_W[-c(1:28),]
                    , control = list(trace = trace, maxit = 150)
                    , method = "SANN")
exp(unname(optim.vals$par))

## We feed the last parameters of SANN in as the first values of Nelder-Mead
optim.vals <- optim(par = optim.vals$par
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = FLUepiCurve_W[-c(1:28),]
                    , control = list(trace = trace, maxit = 800, reltol = 10^-7)
                    , method = "Nelder-Mead"
                    , hessian = T)
optim.vals
MLEfits <- optim.vals$par
MLEfits

exp(unname(MLEfits))

log_alpha.fit <- MLEfits["log_alpha"]
log_alpha.fit
log_Beta.fit <- MLEfits["log_Beta"]
log_Beta.fit
## Look at the output of optim. Understand what it means. Did the algorithm
## converge? Look at ?optim to understand it.

## Plot MLE fit time series
par(bty='n', lwd = 2, las = 1)
fitDat <- simEpidemic(init_state, parms = subsParms(optim.vals$par,  disease_params()))
fitDat
ggplot(fitDat, aes(, x= time, y = daily_deaths)) + 
  geom_line() + 
  geom_line(data = FLUepiCurve_W, aes(x = time, y = ma_7day), col= "green")
with(fitDat, plot(time, daily_deaths, col='blue',lty = ))
points(myDat$time, myDat$sampPrev, col = 'red', pch = 16, cex = 2)
arrows(myDat$time, myDat$uci, myDat$time, myDat$lci, col = 'red', len = .025, angle = 90, code = 3)
legend("topleft", c('truth', 'observed', 'fitted'), lty = c(1, NA, 1), pch = c(NA,16, NA),
       col = c('red', 'red', 'blue'), bty = 'n')
################################################################################











################################################################################
init_state <- c(S = 167304, I = 1, R = 0, D = 0)
simEpidemic <- function(init = init_state, tseq = times, modFunction=SIRmod, parms = disease_params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  simDat$daily_deaths = ceiling(c(0, diff(simDat$D)))
  return(simDat)
}

## Return the negative log of likelihood by convention
nllikelihood <- function(parms = disease_params(), obsDat=FLUepicurve_O[-c(1:28),]) {
  simDat <- simEpidemic(init_state, parms=parms)
  matchedTimes <- simDat$time %in% obsDat$time
  nlls <- -dpois(obsDat$ma_7day, lambda = simDat$daily_deaths[matchedTimes], log = T)
  return(sum(nlls))
}
nllikelihood(disease_params(beta = 3, mu = 0.01))  ## vs some random guessed parameters
#################################################################################
subsParms <- function(fit.params, fixed.params=disease_params()) {
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })
}
## Make likelihood a function of fixed and fitted parameters.
objFXN <- function(fit.params ## paramters to fit
                   , fixed.params =disease_params() ## fixed paramters
                   , obsDat=FLUepicurve_O[-c(1:28),]) {
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat) ## then call likelihood
}
objFXN(guess.params, disease_params())

## Select initial values for fitted parameters from which optimization routine
## will start. If you select bad initial values the algorithm can get stuck on a
## bad set of parameters. You can always try the true values as a starting point
## for this problem, although that's rarely possible in real problems.

init.pars <- c(log_beta = log(0.6), log_mu = log(1/60))
## We will start with SANN optimization since it is stochastic and therefore
## less likely to get stuck in a local minima. But then finish with Nelder-Mead
## optimization which is much faster.

###  NOTE: for trace >0 you see more progress report, bigger numbers show more
###  update
trace <- 3
trace
## SANN: This is stochastic, be CAREFUL -- sometimes it gets stuck at local minima
## for unreasonable parameters. If you see this happen, run it again!
optim.vals <- optim(par = init.pars
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = FLUepicurve_O[-c(1:28),]
                    , control = list(trace = trace, maxit = 150)
                    , method = "SANN")
exp(unname(optim.vals$par))

## We feed the last parameters of SANN in as the first values of Nelder-Mead
optim.vals <- optim(par = optim.vals$par
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = FLUepicurve_O[-c(1:28),]
                    , control = list(trace = trace, maxit = 800, reltol = 10^-7)
                    , method = "Nelder-Mead"
                    , hessian = T)
optim.vals
MLEfits <- optim.vals$par
MLEfits

exp(unname(MLEfits))

log_alpha.fit <- MLEfits["log_alpha"]
log_alpha.fit
log_Beta.fit <- MLEfits["log_Beta"]
log_Beta.fit
## Look at the output of optim. Understand what it means. Did the algorithm
## converge? Look at ?optim to understand it.

## Plot MLE fit time series
par(bty='n', lwd = 2, las = 1)
fitDat <- simEpidemic(init_state, parms = subsParms(optim.vals$par,  disease_params()))
fitDat
ggplot(fitDat, aes(, x= time, y = daily_deaths)) + 
  geom_line() + 
  geom_line(data = FLUepicurve_O, aes(x = time, y = ma_7day), col= "blue")
with(fitDat, plot(time, daily_deaths, col='blue',lty = ))
points(myDat$time, myDat$sampPrev, col = 'red', pch = 16, cex = 2)
arrows(myDat$time, myDat$uci, myDat$time, myDat$lci, col = 'red', len = .025, angle = 90, code = 3)
legend("topleft", c('truth', 'observed', 'fitted'), lty = c(1, NA, 1), pch = c(NA,16, NA),
       col = c('red', 'red', 'blue'), bty = 'n')