#' Update outcome data
#'
#' See comments in gs.
#'
#' @export

outcome_update <- function() {

  dontpanic::msg("Scrape gs")

  outcome_key <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1QyqWuUprTmZObbo2x46G2YrMAxgk2xJENMjRlqp60VQ/edit#gid=0",
                            sheet = "outcome") %>%
    janitor::clean_names()

  dontpanic::msg("Update data")

  usethis::use_data(outcome_key, overwrite = TRUE)

  devtools::load_all()

}
