library(deSolve)
library(tidyverse)

#Initial state
init_state <- c(
  S_w = 63295, I_w = 3, R_w = 0, D_w = 0,
  S_o = 167303, I_o = 1, R_o = 0, D_o = 0
)

#datasets
FittingData_O <- FLUepiCurve_O[-c(1:15),]
FittingData_W <- FLUepiCurve_W[-c(1:15),]

# Times (same as your observed data time range)
times <- seq(0, 160, by = 1)

########################################################################################
# Disease parameters
disease_params <- function(
    beta_oo = 3,
    beta_ow = 0.01,
    beta_ww = 1,
    beta_wo = 0.01,
   # gamma_o = 1 / 8,
   # gamma_w = 1 / 8,
    mu_o = 1/30,
    mu_w = 1/60
   )
  {
  return(as.list(environment()))
}

#######################################################################################

#Creating a function for the SIR model
SIRmod <- function(tt, yy, parms) with(c(parms,as.list(yy)), {
  ## State variables are: S, I, and R
  N_o <- S_o + I_o + R_o
  N_w <- S_w + I_w + R_w
  
  gamma_w = 1/8## total population
  gamma_o = 1/8
  
  ## state variable derivatives (ODE system)
  deriv <- rep(NA,8)
  deriv[1] <- -beta_ww * S_w * I_w / N_w - beta_wo * S_w * I_o/N_o ## Instantaneous rate of change: Susceptibles
  deriv[2] <-	beta_ww * S_w * I_w / N_w + beta_wo * S_w * I_o/N_o - gamma_w * I_w  - mu_w * I_w ## Instantaneous rate of change: Infection class I1
  deriv[3] <- gamma_w * I_w  ##Cumulative recovered people
  deriv[4] <- mu_w * I_w  # Cumulative deaths
  
  deriv[5] <- -beta_oo * S_o * I_o / N_o - beta_ow * S_o * I_w/N_w
  deriv[6] <-  beta_oo * S_o * I_o / N_o + beta_ow * S_o * I_w/N_w- gamma_o * I_o  - mu_o * I_o
  deriv[7] <- gamma_o * I_o
  deriv[8] <- mu_o * I_o # Cumulative deaths
  
  return(list(deriv))
})


#######################################################################################
#simEpidemic() function
simEpidemic <- function(init = init_state, tseq = times, modFunction=SIRmod, parms = disease_params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  simDat$daily_deaths_o = ceiling(c(0, diff(simDat$D_o)))
  simDat$daily_deaths_w = ceiling(c(0, diff(simDat$D_w)))
  return(simDat)
}

#######################################################################################
#nllikelihood() function
nllikelihood <- function(parms = disease_params(), obsDat1 = FittingData_O, obsDat2 = FittingData_W ) {
  simDat <- simEpidemic(init_state, parms = parms)
  #simDat <- simDat %>% filter(!is.na(daily_deaths_total))
  
  matched1 <- simDat$time %in% obsDat1$time
  matched2 <- simDat$time %in% obsDat2$time
  
  if (sum(matched1) == 0||sum(matched2) == 0) stop("No matching time points.")
  
  nlls1 <- -dpois(obsDat1$ma_7day, lambda = simDat$daily_deaths_o[matched1], log = T)
  nlls2 <- -dpois(obsDat2$ma_7day, lambda = simDat$daily_deaths_w[matched2], log = T)
  return(sum(nlls1)+sum(nlls2))
  
}
nllikelihood()

########################################################################################
# Substitute parameters
subsParms <- function(fit.params, fixed.params = disease_params()) {
  #fixed.params <- as.list(fixed.params)  # <-- ADD THIS LINE
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_', '', nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })
}
# subsParms <- function(fit.params, fixed.params=disease_params()) {
#   within(fixed.params, {
#     loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
#     unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
#     for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
#     for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
#     rm(nm, loggedParms, unloggedParms)
#   })
# }
#########################################################################################
# Objective function
objFXN <- function(fit.params ## paramters to fit
                   , fixed.params =disease_params() ## fixed paramters
                   , obsDat1=FittingData_O, obsDat2=FittingData_W) {
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat1 = obsDat1, obsDat2 = obsDat2) ## then call likelihood
}
# guess.params <- c(log_beta = log(0.6), log_mu = log(1/60))
# guess.params
# subsParms(guess.params, disease_params())
# objFXN(guess.params, disease_params())

##########################################################################################
# Initial guess
guess.params <- c(log_beta_oo = log(3),log_beta_ow = log(0.01),log_beta_ww = log(1),log_beta_wo = log(0.01),
                  log_mu_o = log(1/30),log_mu_w = log(1/60) )
init.pars <- guess.params


##########################################################################################
# Stochastic optimization (SANN)
trace <- 3
optim.vals <- optim(
  par = init.pars,
  fn = objFXN,
  fixed.params = disease_params(),
  obsDat1 = FittingData_O, obsDat2=FittingData_W,
  control = list(trace = trace, maxit = 150),
  method = "SANN"
)
##########################################################################################
# Deterministic optimization (Nelder-Mead)
optim.vals <- optim(
  par = optim.vals$par,
  fn = objFXN,
  fixed.params = disease_params(),
  obsDat1 = FittingData_O, obsDat2=FittingData_W,
  control = list(trace = trace, maxit = 800, reltol = 1e-7),
  method = "Nelder-Mead",
  hessian = TRUE
)
########################################################################################
# Maximum Likelihood Estimates
MLEfits <- optim.vals$par
exp(unname(MLEfits))

#FITSWELL <- MLEfits
########################################################################################
# Simulate fitted model and plot
fitDat <- simEpidemic(init_state, parms = subsParms(MLEfits, disease_params()))

#########################################################################################
# Plot
#library(ggplot2)
combinedata <- data.frame(time = c(fitDat$time, FLUepiCurve_O$time), 
                          daily_deaths=c(fitDat$daily_deaths_o, FLUepiCurve_O$ma_7day),
                          output=c(rep("model", nrow(fitDat)), rep("data", nrow(FLUepiCurve_O))))
combinedata <- combinedata %>% left_join(time.data, by="time")

ggplot(combinedata, aes(x = death_date, y = daily_deaths, col=output)) + 
  geom_line(size = 1) + 
  labs(title = "Fitted vs Observed Deaths - Other", y = "Daily Deaths", x = "Time (in days)") +
 theme_minimal()+
  theme(
    plot.title = element_text(size = 18, face = "bold"),     # Title font size
    axis.title.x = element_text(size = 16),                  # X-axis label font size
    axis.title.y = element_text(size = 16),                  # Y-axis label font size
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = c(0.8, 0.8)
  )



combinedata <- data.frame(time = c(fitDat$time, FLUepiCurve_W$time), 
                          daily_deaths=c(fitDat$daily_deaths_w, FLUepiCurve_W$ma_7day),
                          output=c(rep("model", nrow(fitDat)), rep("data", nrow(FLUepiCurve_W))))
combinedata <- combinedata %>% left_join(time.data, by="time")

ggplot(combinedata, aes(x = death_date, y = daily_deaths, col=output)) + 
  geom_line(size = 1) + 
  labs(title = "Fitted vs Observed Deaths - White", y = "Daily Deaths", x = "Time (in days)", lwd = 3) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 18, face = "bold"),     # Title font size
    axis.title.x = element_text(size = 16),                  # X-axis label font size
    axis.title.y = element_text(size = 16),                  # Y-axis label font size
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = c(0.8, 0.8)
  )
 


