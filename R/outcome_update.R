#' Update outcome data
#'
#' See comments in gs.
#'
#' @export

outcome_update <- function() {

  message("Scrape gs")

  # outcome_key <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1QyqWuUprTmZObbo2x46G2YrMAxgk2xJENMjRlqp60VQ/edit#gid=0",
  #                           sheet = "outcome") %>%
  #   janitor::clean_names()

  outcome_url <- "https://docs.google.com/spreadsheets/d/1QyqWuUprTmZObbo2x46G2YrMAxgk2xJENMjRlqp60VQ/edit#gid=0"

  message("Update data")

  # usethis::use_data(outcome_key, overwrite = TRUE)

  count_obs <-
    googlesheets4::read_sheet(outcome_url, "Sperm count") %>%
    janitor::clean_names()

  volume_obs <-
  googlesheets4::read_sheet(outcome_url, "Semen volume") %>%
    janitor::clean_names()

  motility_obs <-
    googlesheets4::read_sheet(outcome_url, "Sperm motility") %>%
    janitor::clean_names()

  morphology_obs <-
    googlesheets4::read_sheet(outcome_url, "Sperm morphology") %>%
    janitor::clean_names()

  usethis::use_data(count_obs, overwrite = TRUE)
  usethis::use_data(motility_obs, overwrite = TRUE)
  usethis::use_data(morphology_obs, overwrite = TRUE)
  usethis::use_data(volume_obs, overwrite = TRUE)

  devtools::build()

}
