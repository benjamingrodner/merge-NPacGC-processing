#Aim, load all sample_metadata.csv file in r, rbind, and check if SampleIDs are unique
library(dplyr)
library(lubridate)

setwd("/Users/sachacoesel/Documents/NPac_v2/NPAc_v2_metadata/sample_metadata")

library(dplyr)
library(purrr)
library(readr)

# Set path to directory with CSVs
path <- "."

# List all csv files
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)

meta <- files %>%
  map_dfr(~ read_csv(.x, show_col_types = FALSE, col_types = cols(.default = "c")))

# Convert numeric columns back where needed
meta <- meta %>%
  mutate(across(c(sequencingID, Station, Cast, Latitude.dec, Longitude.dec,
                  Depth.m, Filter.um, Volume.L), as.numeric))

# Check if SampleID is unique
is_unique <- n_distinct(meta$SampleID) == nrow(meta)

if (is_unique) {
  message("✅ All SampleID values are unique.")
} else {
  message("⚠️ SampleID contains duplicates.")
  
  dupes <- meta %>%
    group_by(SampleID) %>%
    filter(n() > 1) %>%
    arrange(SampleID)
  
  print(dupes)
}
