#' Title
#'
#' @param mod
#' @param m_type
#' @param dir_impr
#'
#' @importFrom ggh4x guide_axis_manual
#' @return
#' @export
#'
#' @examples
#'
#'
nma_forest <- function(mod,
                       m_type = c("lor", "smd"),
                       dir_impr = "upper") {
  dirty_xmas_pal <- list(red = "#b12a1b",
                         green = "#67852e")

  active_pal <- dirty_xmas_pal

  active_pal <-
    if (dir_impr == "lower")
      c(active_pal[2], active_pal[1])
  else
    active_pal

  mod_dat <-
    mod %>%
    # needs to work on model object, not on results
    pluck("network", "agd_arm") %>%
    rename(study = .study, arm = .trt)

  int_n <-
    mod_dat %>%
    mutate(across(c(study, arm), as.character)) %>%
    group_by(arm) %>%
    summarise(int_n = sum(n)) %>%
    arrange(desc(int_n))

  stan_dat <-
    mod %>%
    summary() %>%
    pluck("summary") %>%
    janitor::clean_names() %>%
    dplyr::filter(str_detect(parameter, "^d\\[|^tau")) %>%
    mutate(parameter = str_remove(parameter, "^d\\[") %>%
             str_remove("\\]")) %>%
    rename(arm = parameter,
           ci_lb = "x2_5_percent",
           ci_ub = "x97_5_percent") %>%
    select(arm, mean, ci_lb, ci_ub) %>%
    mutate(
      mod_type = m_type,
      mean = if_else(mod_type == "lor", exp(mean), mean),
      ci_lb = if_else(mod_type == "lor", exp(ci_lb), ci_lb),
      ci_ub = if_else(mod_type == "lor", exp(ci_ub), ci_ub)
    )

  this_tau <- stan_dat %>%
    dplyr::filter(arm == "tau") %>%
    pluck("mean") %>%
    as.numeric() %>%
    round(2) %>%
    as.character()

  n_studies <- mod_dat$study %>% unique() %>% length()

  this_cap <-
    glue(
      "Right y-axis: intervention point estimate [95% credible interval].
      Confidence intervals containing 0 are in grey, significant results in black.
      Point estimates in the direction of improvement, higher, are in green,\n point estimates below no effect in red.\n

Between-intervention heterogeneity: {this_tau}.")

  dir_lgl <- dir_impr == "lower"


  plot_dat <-
    stan_dat %>%
    dplyr::filter(arm != "tau") %>%
    left_join(int_n, by = "arm") %>%
    mutate(
      arm = fct_reorder(arm, mean, .desc = dir_lgl),
      ci = glue::glue("{round(mean, 2)} [{round(ci_lb, 2)}, {round(ci_ub, 2)}]") %>%
        as.character()
      # ,
      # ci = fct_reorder(ci, arm)
    )

  y_axis_first <- levels(plot_dat$arm)
  y_axis_second <- plot_dat %>% arrange(arm) %>% pluck("ci")

  int_ci <- function(int) {
    plot_dat %>%
      dplyr::filter(arm == int) %>%
      pull(ci)
  }

  effect_sizes <-
    tibble(
      effect_size = c(0, 0.2, -0.2, 0.5, -0.5, 0.8, -0.8),
      effect_label = c("no effect", rep("small (0.2)", 2), rep("moderate (0.5)", 2), rep("large (0.8)", 2)) %>%
        fct_relevel("no effect", "small (0.2)",
                    "moderate (0.5)", "large (0.8)")
    ) %>%
    filter(effect_size < max(plot_dat$ci_ub) &
             effect_size > min(plot_dat$ci_lb))

  plot_dat %>%
    ggplot() +
    theme_minimal(base_family = "serif") +
    # geom_vline(
    #   data = effect_sizes,
    #   aes(
    #     xintercept = effect_size,
    #     linetype = effect_label
    #   ),
    #   alpha = 0.3,
    #   size = 0.5
    # ) +
    geom_vline(xintercept = 0,
               alpha = 0.5,
               linetype = "dotted") +
    geom_segment(
      aes(
        x = ci_lb,
        xend = ci_ub,
        y = arm,
        yend = arm,
        colour = I(if_else(ci_lb < 0 &
                             ci_ub > 0, "grey", "black")),
      ),
      alpha = 0.95,
      size = 1.2
    ) +
    geom_point(
      aes(
        x = mean,
        y = arm,
        colour = I(if_else(mean < 0, active_pal[[1]], active_pal[[2]])),
        size = int_n,
        NULL
      ),
      alpha = 0.85,
      # set to diamond
      shape = 18 # http://www.sthda.com/english/wiki/ggplot2-point-shapes
    ) +
    theme(
      # panel.grid = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right",
      legend.direction = "vertical"
    ) +
    guides(y.sec = ggh4x::guide_axis_manual(breaks = y_axis_first,
                                            labels = y_axis_second))  +
    labs(
      x = "Mean difference",
      size = "Intervention sample size",
      y = "",
      caption = this_cap,
      # linetype = "Effect size",
      subtitle = glue("Network meta-analysis on {n_studies} studies with {sum({mod_dat$n})} total participants") %>% str_wrap()
    )



}
