#' Calculate smd accounting for arms
#'
#' Using David's code.
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples

smd_calc <- function(df){
  df %>%
    arrange(outcome, study, desc(control)) %>%
    group_by(outcome, study) %>%
    mutate(
      se = sd/sqrt(n),
      mi = sum(n) - n(),
      cmi = exp(lgamma(mi / 2) - log(sqrt(mi / 2)) - lgamma((mi - 1) /
                                                              2)),
      sdpool = sqrt(weighted.mean(sd ^ 2, n - 1)),
      smd = if_else(control == TRUE, NA_real_, (mean - first(mean)) / sdpool * cmi),
      se_smd = if_else(control == TRUE,
                       se / sdpool * cmi,
                       sqrt((n + first(
                         n
                       )) / (n * first(
                         n
                       )) + smd ^ 2 / (2 * (
                         n + first(n)
                       ))))
    ) %>%
    select(study, intervention, control, smd, se_smd, n, everything())
}
