# Load required library
library(deSolve)
library(dplyr)

# Define the SIR model with mortality
sir_with_mortality <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    
    # Differential equations
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I - mu * I
    dR <- gamma * I
    dD <- mu * I  # Cumulative deaths
    
    list(c(dS, dI, dR, dD))
  })
}

# Initial values
init_state <- c(S = 230599, I = 1, R = 0, D = 0)

# Parameters
params <- c(
  beta = 0.6,   # Infection rate
  gamma = 0.2,  # Recovery rate
  mu = 0.01     # Mortality rate
)

# Time steps
times <- seq(0, 160, by = 1)

# Solve ODE
output <- as.data.frame(ode(y = init_state, times = times, func = sir_with_mortality, parms = params))

# Add daily deaths column (difference between successive cumulative deaths)
output <- output %>%
  mutate(daily_deaths = c(NA, diff(D)))

# Print a preview
head(output, 10)
View(output)
