library(odin2)

span_flu <- odin2::odin({
	
	initial(N) <- N0
	initial(S) <- S0
	# initial(E) <- E0
	initial(I) <- I0
	initial(R) <- R0
	
	initial(new_cases, zero_every = 7) <- 0
	# initial(new_vaccination, zero_every = 7) <- 0
	initial(all_deaths, zero_every = 7) <- 0
  initial(all_recovery, zero_every = 7) <- 0

	## process flows
	#birth <- parameter(10000/365.25)
	beta <- parameter(0.4)
	lam <- beta * I / N
	
	gam <- parameter(1/10) # recovery
	mu_d <- parameter(1/3) # disease death
	
	N0 <- parameter(230600)
	I0 <- 1
	R0 <- 0
	S0 <- N0 - I0 - R0
	
	p_Sout <- 1 - exp(-(lam) * dt)
	p_Ideath <- 1 - exp(-(mu_d)*dt)
	p_Recovery <- 1- exp(-(gam)*dt)
	
	n_Sout <- Binomial(S, p_Sout)
	n_Iout <- Binomial(I, p_Ideath + p_Recovery)
	
	#n_deathsS <- Binomial(n_Sout, mu_b/(lam + tau + mu_b))
	n_deathsI <- Binomial(n_Iout, p_Ideath/(p_Ideath + p_Recovery))
  n_Recover <- n_Iout - n_deathsI

	update(N) <- N - n_deathsI
	update(S) <- S - n_Sout
	update(I) <- I - n_deathsI - n_Recover + n_Sout
	update(R) <- R + n_Recover
	
	update(new_cases) <- new_cases + n_Sout 
	#update(new_vaccination) <- new_vaccination + n_vaccine
	update(all_deaths) <- all_deaths + n_deathsI
	update(all_recovery) <- all_recovery + n_Recover
	
})

span_sys <- dust2::dust_system_create(
	span_flu(),
	pars = list(
	  beta = 0.5,
	  gam = 1 / 10,
	  mu_d = 1 / 3,
	  N0 = 230600,
    dt = 1
	),
	n_particles = 100,
	n_groups = 1,
	seed = NULL,
	deterministic = FALSE,
	n_threads = 8
)

time <- 0:200
dust2::dust_system_set_state_initial(span_sys)

dust2::dust_system_set_time(span_sys, 0)
y <- dust2::dust_system_simulate(span_sys, time )


library(tidyverse)

# Convert to a tidy tibble
y_tidy <- purrr::map_dfr(1:dim(y)[2], function(p) {
  df <- as.data.frame(t(y[, p, ]))
  colnames(df) <- c("N0", "S", "I", "R", "New_cases", "All_deaths", "All_recovery")  # adjust as needed
  df$time <- 0:(nrow(df) - 1)
  df$particle <- p
  df
})
y_tidy

#outbreak_df <- data.frame(t(y))

#names(outbreak_df) <- c("N0", "S", "I", "R" , "New_case" ,"All_deaths", "All_recovery")

#outbreak_df
#create a new column based on the row number, to represent the time step

y_tidy |> ggplot(aes(x = time , y = I, colour = as.factor(particle) )) + geom_line()
y_tidy1 <- y_tidy |> 
  group_by(time) |> 
  summarise(I=mean(I))
view(y_tidy)

y_tidy1 |> ggplot(aes(x = time , y = I) ) + geom_line(col="red")

ggplot(data = y_tidy, aes(x = time , y = I, col = as.factor(particle)  ),) + geom_line(alpha= .3) + 
  geom_line(data = y_tidy1, aes(x = time, y= I), col = "red", size = 1.5) + labs(title = "Infected over time",
                                                                     x = "Time(days)",
                                                                     y = "Number of Infected")

# c(1,2,3,4,5) |> sum()
# sum(c(1,2,3,4,5))

#OUTPUT OF Y
# N
# S
# I
# R
library(dust2)
  ?dust_system_simulate
