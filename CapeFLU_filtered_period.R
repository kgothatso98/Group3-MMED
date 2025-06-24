library(ggplot2)

ggplot(CapeFLU_filtered_period, aes(x = Duration_in_days)) +
  geom_histogram(binwidth = 1, fill = "darkorange", color = "black") +
  scale_x_continuous(limits = c(0, 50)) +
  labs(
    title = "Frequency of Duration of Illness",
    x = "Duration of Illness (Days)",
    y = "Frequency"
  ) +
  theme_minimal()

summary(CapeFLU_filtered_period$Duration_in_days[CapeFLU_filtered_period$RACE=="White"])
summary(CapeFLU_filtered_period$Duration_in_days[CapeFLU_filtered_period$RACE=="Other"])



library(dplyr)
library(stringr)

CapeFLU_filtered <- CapeFLU %>%
  filter(Causes.of.death2=="Pneumonia/Influenza")
#    str_detect(str_to_lower(Causes.of.death), "Spanish Influenza| Pneumonia") &
#      str_detect(str_to_lower(Causes.of.death2), "pneumonia|influenza")
#  )

  
library(dplyr)
library(lubridate)

# Step 1: Create a proper Date column
CapeFLU_filtered <- CapeFLU_filtered %>%
  mutate(
    Full.Date = make_date(Year.of.death.number, Month.of.death.number, Day.of.death.number)
  )

# Step 2: Filter for 1 September 1918 to 31 January 1919
CapeFLU_filtered_period <- CapeFLU_filtered %>%
  filter(Full.Date >= as.Date("1918-09-01") & Full.Date <= as.Date("1919-01-31"))
