library(targets)

#suppressMessages({
library(tidyverse)
library(spunk)
library(janitor)
library(targets)
library(multinma)
library(metafor)
library(glue)
library(gt)
conflicted::conflict_prefer("filter", "dplyr")
options(mc.cores = parallel::detectCores())

# })

# list.files("R", full.names = TRUE) %>% map(source)


list(
  # clean raw data ----------------------------------------------------------
  tar_target(tidy_dat,
             function(dat, this_outcome) {
               dat %>%
                 clean_names() %>%
                 mutate(across(where(is.character), tolower)) %>%
                 mutate(outcome = this_outcome,
                        # get study
                        study = str_extract(study_id, "\\w+\\s\\d+\\w")) %>%
                 rename(class = major_intervention_grouping,
                        intervention = grouped_intervention,
                        moderator = type_of_infertility) %>%
                 select(outcome, study, everything()) %>%
                 select(-starts_with("x"), -intervention_detailed, -study_id)
             }),


  tar_target(
    # Suddenly intervention is being read as a list
    # This problem was first observed 27 of April,
    # so, will likely be the result of a datapoint.
    # We need to move along with the rest of the paper,
    # but it would be really good to have a look into this data
    # discrepency
    count_unlisted,
    count_obs %>%
      filter(!is.na(study_id)) %>%
      mutate(across(
        c(grouped_intervention, major_intervention_grouping),
        as.character
      )) %>%
      mutate(across(
        c(contains("_sd"), contains("_mean"),
          contains("_n")),
        as.numeric
      ))
  ),

  tar_target(wide_obs,
             map2_df(
               list(count_unlisted,
                    volume_obs,
                    motility_obs,
                    morphology_obs),
               list("count", "volume", "motility", "morphology"),
               tidy_dat
             )),


  # class needs to be unique for interventions ---------------------------------


  # identify interventions with more than one class label
  tar_target(
    qa_class,
    wide_obs %>%
      group_by(outcome, intervention) %>%
      summarise(
        n_classes = n_distinct(class),
        classes = unique(class) %>% paste(collapse = "; ")
      ) %>%
      filter(n_classes > 1)
  ),

  # identify most-used class labels
  tar_target(
    int_class,
    wide_obs %>%
      count(intervention, class) %>%
      arrange(intervention, desc(n)) %>%
      group_by(intervention) %>%
      filter(class == first(class)) %>%
      select(-n)
  ),

  # patch in unqiue class labels
  tar_target(
    obs_class_fix,
    wide_obs %>%
      select(-class) %>%
      left_join(int_class) %>%
      select(outcome, intervention, class, moderator, study,
             everything())
  ),


  # set model input ---------------------------------------------------------

  tar_target(model_dat,
             obs_class_fix),

  tar_target(write_model, {
    spunk_dat <- model_dat

    usethis::use_data(spunk_dat, overwrite = TRUE)
  }),



  # attempt multilevel errors (doesn't work) ------------------------------------------------

  tar_target(
    obs_control,
    model_dat %>%
      select(outcome, study, starts_with("control"), class) %>%
      group_by(outcome, study) %>%
      mutate(n_per_study = n_distinct(control_n)) %>%
      sample_n(1) %>%
      select(-n_per_study) %>%
      rename_with( ~ str_remove(.x, "control_")) %>%
      rename(intervention = control) %>%
      mutate(
        control = TRUE,
        class = if_else(str_detect(intervention, "placebo"), "placebo", class),
        class = if_else(
          str_detect(intervention, "supplement"),
          "dietary supplements",
          class
        )
      )

  ),

  tar_target(
    obs_int,
    model_dat %>%
      select(-starts_with("control")) %>%
      rename_with( ~ str_remove(.x, "intervention_")) %>%
      mutate(control = FALSE)

  ),


  tar_target(
    obs_long,
    bind_rows(obs_control, obs_int) %>%
      arrange(outcome, study, desc(control))
  ),

  tar_target(
    check_classes,
    obs_long %>%
      group_by(intervention) %>%
      summarise(
        class_n = n_distinct(class),
        classes = unique(class) %>% paste(collapse = ";")
      ) %>%
      arrange(desc(class_n))
  ),

  tar_target(
    outcome_groups,
    obs_long %>%
      group_by(outcome) %>%
      tar_group(),
    iteration = "group"

  ),

  tar_target(
    fit_dat,
    outcome_groups %>%
      mutate(class = if_else(
        intervention == "vitamin c/e",
        "vitamins",
        class
      )) %>%
      smd_calc() %>%
      mutate(
        class = str_to_sentence(class),
        intervention = str_to_sentence(intervention),
        intervention =
          str_replace(intervention,
                      "Zinc/zinc", "Zinc/Zinc") %>%
          str_replace("d3", "D3") %>%
          str_replace("q10", "Q10") %>%
          str_replace("l-", "L-") %>%
          str_replace("/no", "/No") %>%
          str_replace("c/e", "C/E") %>%
          str_replace("c/", "C/") %>%
          str_replace("ZinC", "Zinc") %>%
          str_replace("combos", "comb") %>%
          str_replace("\\+ala$", "\\+ALA") %>%
          str_replace("\\se", "\\sE") %>%
          str_replace("/min", "/Min"),
        class = str_replace(class, "combos", "comb") %>%
          str_replace("/min", "/Min")
      ),
    pattern = map(outcome_groups),
    iteration = "list"
  ),

  tar_target(
    fit_class_qa,
    fit_dat %>%
      group_by(outcome, intervention) %>%
      summarise(n_class = n_distinct(class)) %>%
      filter(n_class > 1),
    pattern = map(fit_dat),
    iteration = "list"
  ),


  # fit models (using arm-based MD) --------------------------------------------------------------


  tar_target(
    fit_arms,
    set_agd_arm(
      data = fit_dat,
      y = mean,
      se = sd / sqrt(n),
      sample_size = n,
      study = study,
      trt = intervention,
      trt_ref = "Placebo/No treatment",
      trt_class = class

    ) %>%
      nma(trt_effects = "random"),
    pattern = map(fit_dat),
    iteration = "list"
  ),

  # tar_target(
  #   fit_mod,
  #   set_agd_arm(
  #     data = fit_dat,
  #     y = mean,
  #     se = sd/sqrt(n),
  #     sample_size = n,
  #     study = study,
  #     trt = intervention,
  #     trt_ref = "placebo/no treatment",
  #     trt_class = class
  #
  #   ) %>%
  #     nma(trt_effects = "random",
  #         regression = ~.trt:moderator
  #     ),
  #   pattern = map(fit_dat),
  #   iteration = "list"
  # ),

  # tar_target(
  #   fit_nma,
  #   set_agd_contrast(
  #     data = fit_dat,
  #     y = smd,
  #     se = se_smd,
  #     sample_size = n,
  #     study = study,
  #     trt = intervention,
  #     trt_ref = "placebo/no treatment",
  #     trt_class = class
  #
  #   ) %>%
  #     nma(trt_effects = "random",
  #         regression = ~.trt:moderator
  #         ),
  #   pattern = map(fit_dat)
  # ),
  tar_target(
    nma_sucra,
    {
      this_outcome <-
        fit_arms$network$agd_arm$outcome %>% unique()


      multinma::posterior_rank_probs(fit_arms,
                                     lower_better = FALSE,
                                     sucra = TRUE) %>%
        plot() +
        scale_x_continuous(breaks = seq(0, 25, 5)) +
        theme(axis.text = element_text(size = 6),
              strip.text = element_text(size = 5)) +
        labs(
          title = sprintf("SUCRA"),
          subtitle = str_to_sentence(this_outcome),
          x = "Rank"
        )

      ggsave(sprintf("img/%s-sucra.png",
                     this_outcome),
             limitsize = FALSE)
    },
    pattern = map(fit_arms),
    iteration = "list"
  ),


  # extract model estimates -------------------------------------------------

  tar_target(
    nma_rank,
    fit_arms %>%
      posterior_ranks(probs = c(0.025, 0.975),
                      lower_better = FALSE) %>%
      as_tibble() %>%
      clean_names() %>%
      mutate(parameter = str_remove(parameter, "rank\\[") %>%
               str_remove("\\]")) %>%
      select(
        intervention = parameter,
        mean,
        sd,
        lb = x2_5_percent,
        ub = x97_5_percent
      ) %>%
      rename_with(~ glue("{.x}_rank")) %>%
      mutate(outcome = fit_arms$network$agd_arm$outcome %>% unique()) %>%
      select(outcome, everything()),
    pattern = map(fit_arms),
    iteration = "list"
  ),

  tar_target(
    nma_rel,
    fit_arms %>%
      relative_effects(probs = c(0.025, 0.975)) %>%
      as_tibble() %>%
      clean_names() %>%
      mutate(parameter = str_remove(parameter, "d\\[") %>%
               str_remove("\\]")) %>%
      select(
        intervention = parameter,
        mean,
        sd,
        lb = x2_5_percent,
        ub = x97_5_percent
      ) %>%
      rename_with(~ glue("{.x}_rel"))  %>%
      mutate(outcome = fit_arms$network$agd_arm$outcome %>% unique()) %>%
      select(outcome, everything()),
    pattern = map(fit_arms),
    iteration = "list"
  ),

  # wrangle rank-rel matrix -------------------------------------------------
  tar_target(diff_calc,
             function(x, y) {
               diff <- abs(x - y)

               if_else(x < y, diff, -diff)
             }),


  tar_target(
    rank_diff,
    nma_rank %>%
      slice(rep(1:n(), each = n())) %>%
      rename(first_rank = intervention_rank) %>%
      mutate(intervention_rank = nma_rank$intervention_rank %>%
               rep(times = nrow(nma_rank))) %>%
      left_join(nma_rank %>%
                  rename_with(~ glue("{.x}_comp"), -intervention_rank)) %>%
      mutate(
        mean_diff = diff_calc(mean_rank, mean_rank_comp),
        sd_diff = diff_calc(sd_rank, sd_rank_comp),
        lb_diff = diff_calc(lb_rank, lb_rank_comp),
        ub_diff = diff_calc(ub_rank, ub_rank_comp)
      ) %>%
      mutate(
        rank_diff =
          glue("{round(mean_diff)} ({round(lb_diff)}, {round(ub_diff)})")
      ) %>%
      select(
        outcome,
        int_1 = first_rank,
        int_2 = intervention_rank,
        contains("diff"),
        everything()
      ),
    pattern = map(nma_rank),
    iteration = "list"
  ),

  tar_target(
    rel_diff,
    nma_rel %>%
      slice(rep(1:n(), each = n())) %>%
      rename(first_rel = intervention_rel) %>%
      mutate(intervention_rel = nma_rel$intervention_rel %>%
               rep(times = nrow(nma_rel))) %>%
      left_join(nma_rel %>%
                  rename_with(~ glue("{.x}_comp"), -intervention_rel)) %>%
      mutate(
        mean_diff = diff_calc(mean_rel, mean_rel_comp),
        sd_diff = diff_calc(sd_rel, sd_rel_comp),
        lb_diff = diff_calc(lb_rel, lb_rel_comp),
        ub_diff = diff_calc(ub_rel, ub_rel_comp)
      ) %>%
      mutate(
        rel_diff =
          glue(
            "{round(mean_diff, 2)} ({round(lb_diff, 2)}, {round(ub_diff, 2)})"
          )
      ) %>%
      select(
        outcome,
        int_1 = first_rel,
        int_2 = intervention_rel,
        ends_with("diff"),
        everything()
      ),
    pattern = map(nma_rel),
    iteration = "list"
  ),


  # gt rank-rel matrix ------------------------------------------------------

  tar_target(
    rank_rel_dat,
    {
      rels <- rel_diff %>%
        mutate(
          rel_diff = str_replace(rel_diff, "\\s\\(", "<br>\\("),
          rel_diff = if_else(
            int_1 == int_2,
            glue(
              "<b>{round(mean_rel, 2)} <br> ({round(lb_rel, 2)}, {round(ub_rel, 2)})</b>"
            ),
            rel_diff
          )
        ) %>%
        select(outcome, int_1, int_2, rel_diff)

      ranks <- rank_diff %>%
        mutate(
          ub = round(ub_rank),
          rank_diff = if_else(
            int_1 == int_2,
            glue(
              "<b>{round(mean_rank)}  ({round(lb_rank)}, {round(ub_rank)})</b>"
            ),
            rank_diff
          ),
          rank_diff = glue("<i>{rank_diff}</i>")
        ) %>%
        select(outcome, int_1, int_2, rank_diff)


      full_join(rels, ranks) %>%
        mutate(diff = glue("{rel_diff} <br>
      {rank_diff}"))
    },
    pattern = map(rel_diff, rank_diff),
    iteration = "list"
  ),

  tar_target(
    rank_rel_class,
    left_join(rank_rel_dat, int_class,
              by = c("int_1" = "intervention")) %>%
      rename(class_1 = class) %>%
      left_join(int_class,
                by = c("int_2" = "intervention")) %>%
      rename(class_2 = class)
    # mutate(
    #   class_1 = map_chr(int_1, get_int_class),
    #   class_2 = map_chr(int_2, get_int_class)
    # )

    ,
    pattern = map(rank_rel_dat),
    iteration = "list"
  ),

  tar_target(
    rank_rel_class_key,
    rank_rel_class %>%
      select(outcome, class_1) %>%
      distinct() %>%
      rename(class = class_1),
    pattern = map(rank_rel_class),
    iteration = "list"
  ),

  tar_target(
    rank_rel_wide,
    rank_rel_class %>%
      select(-rel_diff, -rank_diff) %>%
      pivot_wider(
        names_from = int_2,
        values_from = diff,
        id_cols = c(outcome, int_1)
      )
    ,
    pattern = map(rank_rel_class),
    iteration = "list"

  ),

  tar_target(rank_rel_gt,
             {
               this_outcome <- rank_rel_wide$outcome %>% unique()

               this_gt <-
                 rank_rel_wide %>%
                 select(-outcome) %>%
                 rename(Intervention = int_1) %>%
                 gt() %>%
                 fmt_markdown(columns = everything()) %>%
                 cols_width(everything() ~ px(60)) %>%
                 tab_style(style = cell_text(weight = "bold"),
                           locations = cells_body(columns = Intervention)) %>%
                 opt_table_lines("all") %>%
                 tab_options(
                   data_row.padding = px(2),
                   data_row.padding.horizontal = px(2),
                   column_labels.font.weight = "bold",
                   table.font.size = 8
                 ) %>%
                 tab_header(
                   title = sprintf(
                     "Difference between %s NMA estimated relative effect compared to placebo and ranks",
                     this_outcome
                   ),
                   subtitle = "Ranks in italics, and bold diagonal is relative effect and rank compared to placebo"
                 )

               gt_path <- sprintf("img/%s-matrix.pdf", this_outcome)

               sprintf("Writing table to:\n%s", gt_path) %>%
                 message()

               write_rds(this_gt,
                         sprintf("img/%s-matrix.rds", this_outcome))
               gtsave(this_gt, gt_path)

             },
             pattern = map(rank_rel_wide)),


  # nma net -----------------------------------------------------------------

  tar_target(plot_network,
             {
               this_outcome <- fit_arms$network$agd_arm$outcome %>% unique()

               fit_arms$network %>%
                 plot(
                   weight_nodes = TRUE,
                   weight_edges = TRUE,
                   show_trt_class = TRUE
                 )  +
                 ggplot2::theme(
                   text = element_text(size = 8),
                   legend.text = element_text(size = 12),
                   legend.position = "bottom",
                   legend.box = "vertical"
                 )

               sprintf("img/%s-net.png", this_outcome) %>%
                 ggsave(
                   filename = .,
                   width = 4200,
                   height = 3000,
                   units = "px"
                 )
             },
             pattern = map(fit_arms)),


  # nma pico ----------------------------------------------------------------


  tar_target(
    pico_dat,
    {
      this_dat <-
        fit_arms$network$agd_arm

      tibble(
        outcome = this_dat$outcome %>% unique(),
        participants = sum(this_dat$n),
        interventions =
          unique(this_dat$.trt) %>%
          paste(collapse = "; ") %>%
          paste0("."),
        classes =
          unique(this_dat$.trtclass) %>%
          paste(collapse = "; ") %>%
          paste0(".")
      )
    },
    pattern = map(fit_arms),
    iteration = "list"
  ),

  tar_target(
    pico_gt,
    {
      this_gt <-
        pico_dat %>%
        t() %>%
        as_tibble(rownames = "cat") %>%
        rename(description = V1) %>%
        gt() %>%
        cols_width(cat ~ px(100)) %>%
        cols_width(description ~ px(300)) %>%
        tab_style(
          style = cell_text(v_align = "top"),
          locations = cells_body(columns = everything(),
                                 rows = everything())
        ) %>%
        tab_style(style = cell_text(weight = "bold", align = "right"),
                  cells_body(columns = cat)) %>%
        opt_table_lines("none") %>%
        tab_options(column_labels.hidden = TRUE)

      img_path <-
        sprintf("img/%s-pico.png", pico_dat$outcome)
      gtsave(this_gt, img_path)
    }
    ,
    pattern = map(pico_dat),
    iteration = "list"
  ),




  # nma sof tab -------------------------------------------------------------

  # paste together the rel
  tar_target(
    sof_rel,
    nma_rel %>%
      arrange(desc(mean_rel)) %>%
      mutate(across(c(
        mean_rel, lb_rel, ub_rel
      ), round, 2)) %>%
      mutate(rel = glue("{mean_rel} ({lb_rel}, {ub_rel})")) %>%
      select(outcome, intervention = intervention_rel, rel),
    pattern = map(nma_rel),
    iteration = "list"

  ),

  # paste rank results
  tar_target(
    sof_rank,
    nma_rank %>%
      mutate(across(c(
        mean_rank, lb_rank, ub_rank
      ), round)) %>%
      mutate(rank = glue("{mean_rank} ({lb_rank}, {ub_rank})")) %>%
      select(outcome, intervention = intervention_rank, rank),
    pattern = map(nma_rank),
    iteration = "list"
  ),

  # get sample sizes
  tar_target(
    sof_n,
    fit_dat %>%
      select(study, intervention, class, n) %>%
      group_by(outcome, intervention, class) %>%
      summarise(n = sum(n)),
    pattern = map(fit_dat),
    iteration = "list"
  ),

  # join rank and rel and sample sizes
  tar_target(
    sof_dat,
    full_join(sof_rel,
              sof_rank) %>%
      full_join(sof_n),
    pattern = map(sof_rel, sof_rank, sof_n),
    iteration = "list"
  ),

  # make gt table
  tar_target(sof_gt, {
    this_outcome <- sof_dat$outcome %>% unique()

    png_path <- sprintf("img/%s-sof.png", this_outcome)

    this_gt <-
      sof_dat %>%
      select(-outcome) %>%
      gt(groupname_col = "class") %>%
      cols_label(
        intervention = "Intervention",
        rel = "Relative effect",
        rank = "Rank",
        n = "Sample size"
      ) %>%
      tab_style(style = cell_text(style = "italic"),
                locations = cells_column_labels(columns = everything())) %>%
      tab_style(style = cell_text(weight = "bold"),
                locations = cells_row_groups()) %>%
      tab_style(style = cell_text(align = "right"),
                locations = cells_body(columns = intervention)) %>%
      tab_style(style = cell_text(align = "left"),
                locations = cells_body(columns = n)) %>%
      tab_style(style = cell_text(align = "right"),
                locations = cells_column_labels(columns = intervention))

    sprintf("Writing image to \n %s", png_path) %>%
      message()

    gtsave(this_gt, png_path)
  },
  pattern = map(sof_dat)),


  # find hth comparisons ----------------------------------------------------

  tar_target(
    hth_comp,
    model_dat %>%
      group_by(outcome, intervention, control) %>%
      summarise(n_studies = n_distinct(study)) %>%
      filter(n_studies > 3)
  ),

  tar_target(
    hth_dat,
    left_join(hth_comp, model_dat),
    pattern = map(hth_comp),
    iteration = "list"
  ),

  tar_target(
    hth_ma,
    escalc(
      measure = "MD",
      m1i = intervention_mean,
      sd1i = intervention_sd,
      n1i = intervention_n,
      m2i = control_mean,
      sd2i = control_sd,
      n2i = control_n,
      data = hth_dat,
      slab = study
    ) %>%
      rma(
        yi = yi,
        vi = vi,
        data = .,
        measure = "MD",
        slab = study
      ),
    pattern = map(hth_dat),
    iteration = "list"
  ),

  # extract results
  tar_target(
    hth_extract_ma,
    tibble(
      mean_ma = hth_ma$b[[1]],
      lb_ma = hth_ma$ci.lb,
      ub_ma = hth_ma$ci.ub
    ),
    pattern = map(hth_ma)
  ),



  # ma remove one -----------------------------------------------------------

  tar_target(
    hth_less,
    {
      study_to_remove <- sample(hth_dat$study, 1)
      hth_dat %>%
        filter(study != study_to_remove)
    },
    pattern = map(hth_dat),
    iteration = "list"
  ),

  tar_target(
    hth_less_ma,
    escalc(
      measure = "MD",
      m1i = intervention_mean,
      sd1i = intervention_sd,
      n1i = intervention_n,
      m2i = control_mean,
      sd2i = control_sd,
      n2i = control_n,
      data = hth_less,
      slab = study
    ) %>%
      rma(
        yi = yi,
        vi = vi,
        data = .,
        measure = "MD",
        slab = study
      ),
    pattern = map(hth_less),
    iteration = "list"
  ),

  # extract results
  tar_target(
    hth_extract_less,
    tibble(
      mean_less = hth_less_ma$b[[1]],
      lb_less = hth_less_ma$ci.lb,
      ub_less = hth_less_ma$ci.ub
    ),
    pattern = map(hth_less_ma)
  ),


  # all ma results ----------------------------------------------------------


  tar_target(
    ma_results_comp,
    bind_cols(hth_comp,
              hth_extract_ma) %>%
      rename_with(~ str_remove(.x, "_ma")) %>%
      mutate(analysis = "MA")
  )
  ,

  tar_target(
    ma_results_less,
    bind_cols(hth_comp, hth_extract_less) %>%
      rename_with(~ str_remove(.x, "_less")) %>%
      mutate(analysis = "MA -1 study")
  ),

  tar_target(
    ma_nma,
    inner_join(
      hth_comp,
      nma_rel %>%
        select(outcome,
               intervention = intervention_rel,
               everything()) %>%
        rename_with( ~ str_remove(.x, "_rel")) %>%
        select(-sd) %>%
        mutate(analysis = "NMA")
    ),
    pattern = map(nma_rel)
  ),

  tar_target(ma_dat,
             bind_rows(ma_results_comp,
                       ma_results_less,
                       ma_nma)),

  tar_target(ma_plot, {
    ma_dat %>%
      ggplot(aes(colour = analysis)) +
      geom_segment(aes(
        x = lb,
        xend = ub,
        y = analysis,
        yend = analysis
      )) +
      geom_point(aes(x = mean, y = analysis)) +
      facet_grid(intervention ~ outcome, scales = "free") +
      theme_minimal() +
      labs(
        title = "Sensitivity analysis: mean point estimates and 95% confidence intervals",
        subtitle = "Comparison of NMA, MA, and MA with one study removed",
        x = "Outcome measure",
        y = "Intervention by MA, MA less one study, and NMA",
        caption = "Analyses were performed on interventions reported by at least 3 studies."
      ) +
      theme(legend.position = "none")

    ggsave("img/sensitivity.png")
  }),

  NULL
)
