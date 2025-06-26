library(tidyverse)
CapeFLU <- readRDS("CapeFLU.rds")
CapeFLU$death_date <- as.Date(CapeFLU$death_date)
CapeFLU$death_date
library(dplyr)

time.data <- data.frame(death_date= seq(as.Date("1918-09-01"),as.Date("1919-01-31"), 1))
time.data$time=0:(nrow(time.data)-1)

FLUepicurve <- CapeFLU %>%
  filter(death_date >= as.Date("1918-09-01") & death_date <= as.Date("1919-01-31")) %>%
  filter(Causes.of.death2 == "Pneumonia/Influenza") %>%
  group_by(death_date) %>%
  summarise(n = n())%>%
  mutate(ma_7day = zoo::rollmean(n, k = 7, fill = NA, align = "center")) %>%
  ungroup() %>% 
  left_join(time.data, by="death_date")
FLUepicurve$ma_7day <- ceiling(FLUepicurve$ma_7day)
FLUepicurve$ma_7day[is.na(FLUepicurve$ma_7day)] <- FLUepicurve$n[is.na(FLUepicurve$ma_7day)]
##################################################################################
#EpiCurve for Race = White
FLUepiCurve_W <- CapeFLU %>%
  filter(death_date >= as.Date("1918-09-01") & death_date <= as.Date("1919-01-31")) %>%
  filter(Causes.of.death2 == "Pneumonia/Influenza") %>%
  filter(RACE == "White") %>%
  group_by(death_date) %>%
  summarise(n = n()) %>%
  mutate(ma_7day = zoo::rollmean(n, k = 7, fill = NA, align = "center")) %>%
  ungroup()%>% 
  left_join(time.data, by="death_date")
FLUepiCurve_W$ma_7day <- ceiling(FLUepiCurve_W$ma_7day)
FLUepiCurve_W$ma_7day[is.na(FLUepiCurve_W$ma_7day)] <- FLUepiCurve_W$n[is.na(FLUepiCurve_W$ma_7day)]
##################################################################################
#EpiCurve for Race = Other
FLUepiCurve_O <- CapeFLU %>%
  filter(death_date >= as.Date("1918-09-01") & death_date <= as.Date("1919-01-31")) %>%
  filter(Causes.of.death2 == "Pneumonia/Influenza") %>%
  filter(RACE == "Other") %>%
  group_by(death_date) %>%
  summarise(n = n()) %>%
  mutate(ma_7day = zoo::rollmean(n, k = 7, fill = NA, align = "center")) %>%
  ungroup()%>% 
  left_join(time.data, by="death_date")
FLUepiCurve_O$ma_7day <- ceiling(FLUepiCurve_O$ma_7day)
FLUepiCurve_O$ma_7day[is.na(FLUepiCurve_O$ma_7day)] <- FLUepiCurve_O$n[is.na(FLUepiCurve_O$ma_7day)]




View(FLUepicurve)
ggplot(FLUepicurve, aes(x = death_date, y = n)) +
  geom_col(fill = "tomato") +
  geom_line(aes(x=death_date, y = ma_7day)) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  labs(
    x = "Months",
    y = "Number of Deaths",
    title = "Daily Pneumonia/Influenza Deaths"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )
FLUepicurve <- FLUepicurve %>%
  mutate(EpiWeek=isoweek(death_date))
FLUepicurve <- FLUepicurve %>%
  # mutate(EpiYear=year(death_date)) %>%
  mutate(
    ma_7day = zoo::rollmean(n, k = 7, fill = NA, align = "center")
  ) %>%
  ungroup()


###################################################################################
FLUepicurve_2 <- CapeFLU %>%
  filter(death_date >= as.Date("1918-09-01") & death_date <= as.Date("1919-01-31")) %>%
  filter(Causes.of.death2 == "Pneumonia/Influenza") %>%
  group_by(RACE, death_date) %>%
  summarise(n = n())


ggplot(FLUepicurve_2, aes(x = death_date, y = n, fill = RACE)) +
  geom_col() +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  labs(
    x = "Months",
    y = "Number of  Weekly Deaths",
    title = "Daily Pneumonia/Influenza Deaths"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )
#################################################################################
#Changing the word "RACE" to "Race" so we have the same column named with 
#the same name
FLUepicurve_2 <- FLUepicurve_2 %>%
  rename(Race = RACE)
View(FLUepicurve_2)
################################################################################
#combined/merge the 2 datasets such that we can have a single dataset
merged_data <- FLUepicurve_2 %>%
  left_join(data_Race, by = "Race")
merged_data
#divide the number of deaths per day by the population of each group
merged_data <- merged_data %>%
  mutate(Death_prop = (n / Population) * 100000)
merged_data
View(merged_data)

ggplot(merged_data, aes(x = death_date, y = Death_prop, color = Race)) +
  geom_col() +
  labs(
    x = "Date",
    y = "Deaths per 100,000 Population",
    title = "Daily Death Rates by Race Group"
  ) +
  theme_minimal()


ggplot(merged_data, aes(x = Death_prop, fill = Race)) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Deaths per 100,000 Population",
    y = "Density",
    title = "Smoothed Distribution of Daily Death Rates by Race"
  ) +
  scale_fill_manual(values = c("White" = "black", "Other" = "lightblue")) +
  theme_minimal()
#################################################################################
#add population to daily deaths
FLUepicurve_2 <- FLUepicurve_2 %>%
  left_join(data_Race, by="Race")

FLUepicurve_2 <- FLUepicurve_2 %>% 
  mutate(Death_prop = 100000*(n/pop.size))
ggplot(FLUepicurve_2, aes(x = death_date, y = Death_prop, fill = Race)) +
  geom_col() +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  labs(
    x = "Months",
    y = "Deaths Per 100000",
    title = "Daily Pneumonia/Influenza Deaths"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )
#################################################################################