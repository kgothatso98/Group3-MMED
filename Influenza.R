library(tidyverse)
CapeFLU <- readRDS("CapeFLU.rds")
CapeFLU$death_date <- as.Date(CapeFLU$death_date)
CapeFLU$death_date
library(dplyr)
FLUepicurve <- CapeFLU %>%
  filter(death_date >= as.Date("1918-09-01") & death_date <= as.Date("1919-01-31")) %>%
  filter(Causes.of.death2 == "Pneumonia/Influenza") %>%
  group_by(death_date) %>%
  summarise(n = n())
View(FLUepicurve)
ggplot(FLUepicurve, aes(x = death_date, y = n)) +
  geom_col(fill = "tomato") +
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
  mutate(EpiYear=year(death_date))
FluepiWeek <- FLUepicurve %>% 
  group_by(EpiYear, EpiWeek) %>% 
  summarise(n = sum(n)) %>% 
  ungroup()
FluepiWeek <- FluepiWeek %>% mutate(Weeknum = seq(1,24,1))
View(FluepiWeek)


FluepiWeek$Year <- ifelse(FluepiWeek$Weeknum <= 19, "1918", "1919")
FluepiWeek$WeekLabel <- paste(FluepiWeek$Year, "- W", FluepiWeek$Weeknum, sep = "")
FluepiWeek$WeekLabel <- factor(FluepiWeek$WeekLabel, levels = FluepiWeek$WeekLabel)
ggplot(FluepiWeek, aes(x = WeekLabel, y = n, fill = Year)) +
  geom_col() +
  scale_fill_manual(values = c("1918" = "lightblue", "1919" = "blue")) +
  labs(
    x = "Week",
    y = "Number of Weekly Deaths",
    fill = "Year",
    title = "Pneumonia/Influenza Weekly Deaths (1918–1919)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)
  )
