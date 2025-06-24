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