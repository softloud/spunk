library(googlesheets4)
library(usethis)

# set url

sperm_url <- "https://docs.google.com/spreadsheets/d/1QyqWuUprTmZObbo2x46G2YrMAxgk2xJENMjRlqp60VQ/edit#gid=0"

# read dat

count_obs <-
  read_sheet(sperm_url, sheet = "Sperm count")

volume_obs <-
  read_sheet(sperm_url, "Semen volume")

morphology_obs <-
  read_sheet(sperm_url, "Sperm morphology")

motility_obs <-
  read_sheet(sperm_url, "Sperm motility")

# write dat
use_data(count_obs, overwrite = TRUE)
use_data(volume_obs, overwrite = TRUE)
use_data(morphology_obs, overwrite = TRUE)
use_data(motility_obs, overwrite = TRUE)
