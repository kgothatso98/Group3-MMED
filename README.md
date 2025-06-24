Modeling the 1918 Spanish Flu in South Africa This repository presents a work-in-progress project exploring the dynamics of the 1918 Spanish Flu epidemic in South Africa using stochastic epidemiological models. It focuses on estimating the basic reproduction number (R₀), accounting for latent infection periods, natural and disease-related mortality, and fitting compartmental models to historical outbreak data. Materials will include code, datasets, and documentation to support research into epidemic heterogeneity and disease dynamics in the South African context.

https://drive.google.com/drive/folders/1z4ZtiZyZr3-RnfMeqN2Ngg3cMG1Emlip?usp=drive_link
##The Code for negative binomial distribution using D fitting with the other D

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
times <- seq(0, 180, by = 1)

# Solve ODE
output <- as.data.frame(ode(y = init_state, times = times, func = sir_with_mortality, parms = params))

# Add daily deaths column (difference between successive cumulative deaths)
output <- output %>%
  mutate(daily_deaths = c(NA, diff(D, lag =1)))
output$daily_deaths[output$time==0]=0
# Print a preview
head(output, 10)
View(output)
#############################################################################################
# Define the SIR model with mortality
sir_with_mortality <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I - mu * I
    dR <- gamma * I
    dD <- mu * I
    
    list(c(dS, dI, dR, dD))
  })
}

# Initial values and parameters
init_state <- c(S = 230599, I = 1, R = 0, D = 0)
params <- c(beta = 0.6, gamma = 0.2, mu = 0.01)
times <- seq(0, 180, by = 1)

# Solve the ODE system
output <- as.data.frame(ode(y = init_state, times = times, func = sir_with_mortality, parms = params))

# Add daily deaths column
output <- output %>%
  mutate(daily_deaths = c(0, diff(D)))

# Group every 30 days
output_grouped <- output %>%
  mutate(period = floor(time / 10)) %>%
  group_by(period) %>%
  summarise(
    period_start = min(time),
    period_end = max(time),
    cum_deaths = max(D)
  ) %>%
  ungroup()

# View grouped summary
print(output_grouped)

# Fit a spline model to the grouped data
spline_fit <- spline(x = output_grouped$period_end, y = output_grouped$cum_deaths, xout = output$time)

# Add the fitted values to the main output
output$D_fitted <- spline_fit$y

# Plot original vs fitted
ggplot(output, aes(x = time)) +
  geom_line(aes(y = D), color = "blue", size = 1, linetype = "solid") +
  geom_line(aes(y = D_fitted), color = "red", size = 1, linetype = "dashed") +
  labs(
    title = "Cumulative Deaths: Original vs 30-Day Grouped Fit",
    y = "Cumulative Deaths",
    x = "Time (days)"
  ) +
  theme_minimal()

######################################################################################

sir_with_mortality <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I - mu * I
    dR <- gamma * I
    dD <- mu * I
    list(c(dS, dI, dR, dD))
  })
}

# --- 2. Initial State, Time, Parameters ---
init_state <- c(S = 230599, I = 1, R = 0, D = 0)
times <- seq(0, 180, by = 1)
true_params <- c(beta = 0.6, gamma = 0.2, mu = 0.01)

# --- 3. Simulate True Model ---
sim_true <- as.data.frame(ode(y = init_state, times = times, func = sir_with_mortality, parms = true_params))
D_true <- sim_true$D
obs_days <- seq(30, 180, by = 30)
true_D_sampled <- D_true[obs_days + 1]  # +1 because R is 1-indexed

# --- 4. Simulate Noisy Observations with Negative Binomial ---
set.seed(123)
dispersion <- 10
observed_D <- rnbinom(length(true_D_sampled), mu = true_D_sampled, size = dispersion)
fit_data <- data.frame(time = obs_days, D_obs = observed_D)

# --- 5. Grid Search for Best Fit ---
beta_vals <- seq(0.5, 0.7, length.out = 5)
gamma_vals <- seq(0.15, 0.25, length.out = 5)
mu_vals <- seq(0.008, 0.015, length.out = 5)

best_sse <- Inf
best_params <- c()
best_fit <- NULL

for (b in beta_vals) {
  for (g in gamma_vals) {
    for (m in mu_vals) {
      params <- c(beta = b, gamma = g, mu = m)
      sim <- as.data.frame(ode(y = init_state, times = times, func = sir_with_mortality, parms = params))
      model_D <- sim$D[obs_days + 1]
      sse <- sum((observed_D - model_D)^2)
      
      if (sse < best_sse) {
        best_sse <- sse
        best_params <- params
        best_fit <- sim
      }
    }
  }
}

# --- 6. Plot Results ---
ggplot() +
  geom_line(data = sim_true, aes(x = time, y = D), color = "red", linetype = "dashed", size = 1) +
  geom_line(data = best_fit, aes(x = time, y = D), color = "blue", size = 1) +
  geom_point(data = fit_data, aes(x = time, y = D_obs), size = 3, color = "black") +
  labs(title = paste0("Best Fit: beta=", round(best_params["beta"],3),
                      ", gamma=", round(best_params["gamma"],3),
                      ", mu=", round(best_params["mu"],3)),
       x = "Time (days)", y = "Cumulative Deaths") +
  theme_minimal()
##############################################################################################
# Define the SIR model with mortality
sir_with_mortality <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I - mu * I
    dR <- gamma * I
    dD <- mu * I
    
    list(c(dS, dI, dR, dD))
  })
}

# Initial values and parameters
init_state <- c(S = 230599, I = 1, R = 0, D = 0)
params <- c(beta = 0.6, gamma = 0.2, mu = 0.01)
times <- seq(0, 180, by = 1)

# Solve the ODE system
output <- as.data.frame(ode(y = init_state, times = times, func = sir_with_mortality, parms = params))

# Add daily deaths column
output <- output %>%
  mutate(daily_deaths = c(0, diff(D)))

# Group every 10 days and summarise cumulative deaths
output_grouped <- output %>%
  mutate(period = floor(time / 10)) %>%
  group_by(period) %>%
  summarise(
    period_start = min(time),
    period_end = max(time),
    cum_deaths = max(D)
  ) %>%
  ungroup()

# View grouped summary
print(output_grouped)

# Create arrow segments from one point to the next
arrows_df <- output_grouped %>%
  mutate(
    xend = lead(period_end),
    yend = lead(cum_deaths)
  ) %>%
  filter(!is.na(xend))  # Remove last row where there's no next point

# Plot using points and arrows
ggplot() +
  geom_point(data = output_grouped, aes(x = period_end, y = cum_deaths), color = "blue", size = 3) +
  geom_segment(data = arrows_df, aes(x = period_end, y = cum_deaths, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm")), color = "red", size = 1) +
  labs(
    title = "Cumulative Deaths by Period with Directional Arrows",
    x = "Time (days)",
    y = "Cumulative Deaths"
  ) +
  theme_minimal()
#####################################################################################
# Define the SIR model with mortality
sir_with_mortality <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I - mu * I
    dR <- gamma * I
    dD <- mu * I
    list(c(dS, dI, dR, dD))
  })
}

# Parameters and initial values
init_state <- c(S = 230599, I = 1, R = 0, D = 0)
params <- c(beta = 0.6, gamma = 0.2, mu = 0.01)
times <- seq(0, 180, by = 1)

# Solve ODE
output <- as.data.frame(ode(y = init_state, times = times, func = sir_with_mortality, parms = params))

# Group by every 10 days
output_grouped <- output %>%
  mutate(period = floor(time / 10)) %>%
  group_by(period) %>%
  summarise(
    period_mid = mean(time),
    cum_deaths = max(D)
  ) %>%
  ungroup()

# Simulate error bars (e.g., ±10% for illustration)
output_grouped <- output_grouped %>%
  mutate(
    lower = cum_deaths * 0.9,
    upper = cum_deaths * 1.1
  )

# Fit a smooth curve (e.g., loess)
smooth_fit <- loess(cum_deaths ~ period_mid, data = output_grouped)
smooth_vals <- predict(smooth_fit, newdata = data.frame(period_mid = output$time))

# Add to original output for plotting
output$smoothed <- smooth_vals

# Plot like the reference image
ggplot() +
  geom_line(data = output, aes(x = time, y = smoothed), color = "red", size = 1) +
  geom_point(data = output_grouped, aes(x = period_mid, y = cum_deaths), color = "red", size = 3) +
  geom_errorbar(data = output_grouped, aes(x = period_mid, ymin = lower, ymax = upper), color = "red", width = 1) +
  labs(
    x = "Time (days)",
    y = "Cumulative Deaths"
  ) +
  theme_minimal()
