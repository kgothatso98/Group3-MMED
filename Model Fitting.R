# Load required library
library(deSolve)
library(dplyr)

# Initial values
Pop1911 <- 167303 #230600 for all; 63296 for W; 167303 for O
init_state <- c(S = Pop1911, I = 1, R = 0, D = 0)
FittingData <- FLUepicurve_O[-c(1:6),]
# Parameters
# params <- c(
#   beta = 0.6,   # Infection rate
#   gamma = 1/8,  # Recovery rate
#   mu = 0.01     # Mortality rate
# )

################################################################################
# Times (same as your observed data time range)
times <- seq(0, 160, by = 1)

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
simEpidemic <- function(init = init_state, tseq = times, modFunction=SIRmod, parms = disease_params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  simDat$daily_deaths = ceiling(c(0, diff(simDat$D)))
  return(simDat)
}

## Return the negative log of likelihood by convention
nllikelihood <- function(parms = disease_params(), obsDat=FittingData) {
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
                   , obsDat=FittingData) {
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
                    , obsDat = FittingData
                    , control = list(trace = trace, maxit = 150)
                    , method = "SANN")
exp(unname(optim.vals$par))

## We feed the last parameters of SANN in as the first values of Nelder-Mead
optim.vals <- optim(par = optim.vals$par
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = FittingData
                    , control = list(trace = trace, maxit = 800, reltol = 10^-7)
                    , method = "Nelder-Mead"
                    , hessian = T)
optim.vals
MLEfits <- optim.vals$par
MLEfits

exp(unname(MLEfits))

## Plot MLE fit time series
par(bty='n', lwd = 2, las = 1)
fitDat <- simEpidemic(init_state, parms = subsParms(optim.vals$par,  disease_params()))
fitDat
ggplot(fitDat, aes(x= time, y = daily_deaths)) + 
  geom_line() + 
  geom_line(data = FittingData, aes(x = time, y = ma_7day), col= "red")
#with(fitDat, plot(time, daily_deaths, col='blue',lty = ))
#points(myDat$time, myDat$sampPrev, col = 'red', pch = 16, cex = 2)
#arrows(myDat$time, myDat$uci, myDat$time, myDat$lci, col = 'red', len = .025, angle = 90, code = 3)
legend("topleft", c('truth', 'observed', 'fitted'), lty = c(1, NA, 1), pch = c(NA,16, NA),
       col = c('red', 'red', 'blue'), bty = 'n')
# ##################################################################################
# # Initial state for the White population
# init_state <- c(S = 63296, I = 1, R = 0, D = 0)
# # Your observed data (must include time and daily_deaths columns)
# # obsDat <- read.csv("your_data.csv")
# # Solve ODE
# # Add daily deaths column (difference between successive cumulative deaths)
# 
# #################################################################################
# #################################################################################
# simEpidemic <- function(init = init_state, tseq = times, modFunction=SIRmod, parms = disease_params()) {
#   simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
#   simDat$daily_deaths = ceiling(c(0, diff(simDat$D)))
#   return(simDat)
# }
# 
# ## Return the negative log of likelihood by convention
# nllikelihood <- function(parms = disease_params(), obsDat=FittingData) {
#   simDat <- simEpidemic(init_state, parms=parms)
#   matchedTimes <- simDat$time %in% obsDat$time
#   nlls <- -dpois(obsDat$ma_7day, lambda = simDat$daily_deaths[matchedTimes], log = T)
#   return(sum(nlls))
#         }
# nllikelihood(disease_params(beta = 3, mu = 0.01))  ## vs some random guessed parameters
# #################################################################################
# subsParms <- function(fit.params, fixed.params=disease_params()) {
#   within(fixed.params, {
#     loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
#     unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
#     for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
#     for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
#     rm(nm, loggedParms, unloggedParms)
#   })
# }
# ## Make likelihood a function of fixed and fitted parameters.
# objFXN <- function(fit.params ## paramters to fit
#                    , fixed.params =disease_params() ## fixed paramters
#                    , obsDat=FittingData) {
#   parms <- subsParms(fit.params, fixed.params)
#   nllikelihood(parms, obsDat = obsDat) ## then call likelihood
# }
# objFXN(guess.params, disease_params())
# 
# ## Select initial values for fitted parameters from which optimization routine
# ## will start. If you select bad initial values the algorithm can get stuck on a
# ## bad set of parameters. You can always try the true values as a starting point
# ## for this problem, although that's rarely possible in real problems.
# 
# init.pars <- c(log_beta = log(0.6), log_mu = log(1/60))
# ## We will start with SANN optimization since it is stochastic and therefore
# ## less likely to get stuck in a local minima. But then finish with Nelder-Mead
# ## optimization which is much faster.
# 
# ###  NOTE: for trace >0 you see more progress report, bigger numbers show more
# ###  update
# trace <- 3
# trace
# ## SANN: This is stochastic, be CAREFUL -- sometimes it gets stuck at local minima
# ## for unreasonable parameters. If you see this happen, run it again!
# optim.vals <- optim(par = init.pars
#                     , objFXN
#                     , fixed.params = disease_params()
#                     , obsDat = FittingData
#                     , control = list(trace = trace, maxit = 150)
#                     , method = "SANN")
# exp(unname(optim.vals$par))
# 
# ## We feed the last parameters of SANN in as the first values of Nelder-Mead
# optim.vals <- optim(par = optim.vals$par
#                     , objFXN
#                     , fixed.params = disease_params()
#                     , obsDat = FittingData
#                     , control = list(trace = trace, maxit = 800, reltol = 10^-7)
#                     , method = "Nelder-Mead"
#                     , hessian = T)
# optim.vals
# MLEfits <- optim.vals$par
# MLEfits
# 
# exp(unname(MLEfits))
# 
# ## Look at the output of optim. Understand what it means. Did the algorithm
# ## converge? Look at ?optim to understand it.
# 
# ## Plot MLE fit time series
# par(bty='n', lwd = 2, las = 1)
# fitDat <- simEpidemic(init_state, parms = subsParms(optim.vals$par,  disease_params()))
# fitDat
# ggplot(fitDat, aes(, x= time, y = daily_deaths)) + 
#   geom_line() + 
#   geom_line(data = FLUepiCurve_W, aes(x = time, y = ma_7day), col= "green")
# with(fitDat, plot(time, daily_deaths, col='blue',lty = ))
# points(myDat$time, myDat$sampPrev, col = 'red', pch = 16, cex = 2)
# arrows(myDat$time, myDat$uci, myDat$time, myDat$lci, col = 'red', len = .025, angle = 90, code = 3)
# legend("topleft", c('truth', 'observed', 'fitted'), lty = c(1, NA, 1), pch = c(NA,16, NA),
#        col = c('red', 'red', 'blue'), bty = 'n')
# ################################################################################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ################################################################################
# init_state <- c(S = 167303, I = 1, R = 0, D = 0)
# simEpidemic <- function(init = init_state, tseq = times, modFunction=SIRmod, parms = disease_params()) {
#   simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
#   simDat$daily_deaths = ceiling(c(0, diff(simDat$D)))
#   return(simDat)
# }
# 
# ## Return the negative log of likelihood by convention
# nllikelihood <- function(parms = disease_params(), obsDat=FLUepicurve_O[-c(1:28),]) {
#   simDat <- simEpidemic(init_state, parms=parms)
#   matchedTimes <- simDat$time %in% obsDat$time
#   nlls <- -dpois(obsDat$ma_7day, lambda = simDat$daily_deaths[matchedTimes], log = T)
#   return(sum(nlls))
# }
# nllikelihood(disease_params(beta = 3, mu = 0.01))  ## vs some random guessed parameters
# #################################################################################
# subsParms <- function(fit.params, fixed.params=disease_params()) {
#   within(fixed.params, {
#     loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
#     unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
#     for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
#     for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
#     rm(nm, loggedParms, unloggedParms)
#   })
# }
# ## Make likelihood a function of fixed and fitted parameters.
# objFXN <- function(fit.params ## paramters to fit
#                    , fixed.params =disease_params() ## fixed paramters
#                    , obsDat=FLUepicurve_O[-c(1:28),]) {
#   parms <- subsParms(fit.params, fixed.params)
#   nllikelihood(parms, obsDat = obsDat) ## then call likelihood
# }
# objFXN(guess.params, disease_params())
# 
# ## Select initial values for fitted parameters from which optimization routine
# ## will start. If you select bad initial values the algorithm can get stuck on a
# ## bad set of parameters. You can always try the true values as a starting point
# ## for this problem, although that's rarely possible in real problems.
# 
# init.pars <- c(log_beta = log(0.6), log_mu = log(1/60))
# ## We will start with SANN optimization since it is stochastic and therefore
# ## less likely to get stuck in a local minima. But then finish with Nelder-Mead
# ## optimization which is much faster.
# 
# ###  NOTE: for trace >0 you see more progress report, bigger numbers show more
# ###  update
# trace <- 3
# trace
# ## SANN: This is stochastic, be CAREFUL -- sometimes it gets stuck at local minima
# ## for unreasonable parameters. If you see this happen, run it again!
# optim.vals <- optim(par = init.pars
#                     , objFXN
#                     , fixed.params = disease_params()
#                     , obsDat = FLUepicurve_O[-c(1:28),]
#                     , control = list(trace = trace, maxit = 150)
#                     , method = "SANN")
# exp(unname(optim.vals$par))
# 
# ## We feed the last parameters of SANN in as the first values of Nelder-Mead
# optim.vals <- optim(par = optim.vals$par
#                     , objFXN
#                     , fixed.params = disease_params()
#                     , obsDat = FLUepicurve_O[-c(1:28),]
#                     , control = list(trace = trace, maxit = 800, reltol = 10^-7)
#                     , method = "Nelder-Mead"
#                     , hessian = T)
# optim.vals
# MLEfits <- optim.vals$par
# MLEfits
# 
# exp(unname(MLEfits))
# 
# log_alpha.fit <- MLEfits["log_alpha"]
# log_alpha.fit
# log_Beta.fit <- MLEfits["log_Beta"]
# log_Beta.fit
# ## Look at the output of optim. Understand what it means. Did the algorithm
# ## converge? Look at ?optim to understand it.
# 
# ## Plot MLE fit time series
# par(bty='n', lwd = 2, las = 1)
# fitDat <- simEpidemic(init_state, parms = subsParms(optim.vals$par,  disease_params()))
# fitDat
# ggplot(fitDat, aes(, x= time, y = daily_deaths)) + 
#   geom_line() + 
#   geom_line(data = FLUepicurve_O, aes(x = time, y = ma_7day), col= "blue")
# with(fitDat, plot(time, daily_deaths, col='blue',lty = ))
# points(myDat$time, myDat$sampPrev, col = 'red', pch = 16, cex = 2)
# arrows(myDat$time, myDat$uci, myDat$time, myDat$lci, col = 'red', len = .025, angle = 90, code = 3)
# legend("topleft", c('truth', 'observed', 'fitted'), lty = c(1, NA, 1), pch = c(NA,16, NA),
#        col = c('red', 'red', 'blue'), bty = 'n')