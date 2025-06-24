library(lubridate)

# Example vector
durations <- c("3 weeks", "8 days", "2 months", "1 year", "one month")

# Use lubridate's duration() on clean numeric values
# But first convert words like "one" to 1 (lubridate can't handle words)
library(stringr)

convert_to_days_lubridate <- function(x) {
  x <- str_to_lower(x)
  x <- str_replace_all(x, "one", "1")
  x <- str_replace_all(x, "two", "2")
  x <- str_replace_all(x, "three", "3")
  x <- str_replace_all(x, "four", "4")
  x <- str_replace_all(x, "five", "5")
  x <- str_replace_all(x, "six", "6")
  x <- str_replace_all(x, "seven", "7")
  x <- str_replace_all(x, "eight", "8")
  x <- str_replace_all(x, "nine", "9")
  
  # Try to parse using lubridate duration
  d <- tryCatch(duration(x), error = function(e) NA)
  
  # Return duration in days
  if (!is.na(d)) {
    return(as.numeric(d, "days"))
  } else {
    return(NA)
  }
}

# Apply to your data column
CapeFLU$Duration_in_days <- sapply(CapeFLU$Duration.of.last.illness, convert_to_days_lubridate)
# Round to nearest whole day
CapeFLU$Duration_in_days <- round(CapeFLU$Duration_in_days)

# Or always round down
# CapeFLU$Duration_in_days <- floor(CapeFLU$Duration_in_days)

# Check result
head(CapeFLU[, c("Duration.of.last.illness", "Duration_in_days")])

