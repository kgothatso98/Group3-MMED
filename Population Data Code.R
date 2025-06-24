library(dplyr)
library(ggplot2)
library(scales)

# Load data
data <- readRDS("Population1911.rds")

# Create consolidated age groups
data <- data %>%
  mutate(
    AgeGroup = case_when(
      Age %in% c("Under 1 Year", "1 Year", "2 Years", "3 Years", "4 Years") ~ "Under 5 Years",
      Age == "5- 9 Years" ~ "5-9 Years",
      Age == "10- 14 Years" ~ "10-14 Years",
      Age == "15- 19 Years" ~ "15-19 Years",
      Age == "20- 24 Years" ~ "20-24 Years",
      Age == "25- 34 Years" ~ "25-34 Years",
      Age == "35- 49 Years" ~ "35-49 Years",
      Age == "50- 69 Years" ~ "50-69 Years",
      Age == "70 Years and Over" ~ "70+ Years",
      TRUE ~ as.character(Age)
    )
  )

# Define proper age group ordering
age_order <- c("Under 5 Years", "5-9 Years", "10-14 Years", "15-19 Years",
               "20-24 Years", "25-34 Years", "35-49 Years", "50-69 Years", "70+ Years")

# Process age data with proper ordering
pop_by_age <- data %>%
  group_by(AgeGroup) %>%
  summarise(pop_size = sum(pop.size, na.rm = TRUE)) %>%
  mutate(AgeGroup = factor(AgeGroup, levels = age_order)) %>%
  arrange(AgeGroup)

# Create enhanced age distribution plot
age_plot <- ggplot(pop_by_age, aes(x = AgeGroup, y = pop_size)) +
  geom_col(fill = "#1f77b4", width = 0.7) +
  geom_text(aes(label = comma(pop_size)), vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = comma, 
                     expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Population Distribution by Age Group\nSouth Africa, 1911 Census",
       x = "Age Group",
       y = "Population Count",
       caption = "Source: South African Census 1911") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(hjust = 1, margin = margin(t = 10))
  )

# Print total population
total_pop <- sum(pop_by_age$pop_size)
cat("Total population in 1911 census:", comma(total_pop), "\n")

# Show the plot
print(age_plot)

# Save the plot (optional)
ggsave("SA_1911_population_by_age_consolidated.png", age_plot, 
       width = 10, height = 6, dpi = 300)




data_Race <- data %>% 
  group_by(Race) %>%
  summarise(pop.size = sum(pop.size))
ggplot(data_Race, aes(x= Race, y = pop.size))+
  geom_col()+
  scale_fill_manual(values = c("Other" = "lightblue", "White"="black"))+
  labs(
    x= "Race",
    y= "Number of Deaths",
    title = "Number of Deaths by Race"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)
  )
##################################################################################

data_Race$Population <- c(167304, 63297)  # Order must match 'Race'
data_Race$death_rate <- (FluepiWeek$n / data_Race$Population) * 100000

library(ggplot2)

ggplot(data_Race, aes(x = Race, y = death_rate, fill = Race)) +
  geom_col() +
  scale_fill_manual(values = c("White" = "black", "Other" = "lightblue")) +
  labs(
    x = "Race Group",
    y = "Deaths per 100,000 Population",
    title = "Death Rate by Race Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 10)
  )