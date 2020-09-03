# Predictive accuracy of various short term death count forecasts
# Jonas Schöley

# Init ------------------------------------------------------------

library(tidyverse)

source('code/00-init.R')

dat <- list()
fig <- list()

# Load fitted models, test and training data ----------------------

load('out/model_fits.RData')

# Residuals by week and sex ---------------------------------------

dat$residuals_by_week_and_sex <-
  mod %>%
  unnest(predictions) %>%
  group_by(cv_id, model_name, sample,
           iso_week, sex) %>%
  summarise(
    observed_deaths =
      sum(observed_deaths, na.rm = T),
    predicted_deaths =
      sum(predicted_deaths, na.rm = T),
    exposure =
      sum(exposure)
  ) %>%
  ungroup() %>%
  mutate(
    observed_mortality =
      observed_deaths/exposure,
    predicted_mortality =
      predicted_deaths/exposure,
    residual_deaths_e =
      observed_deaths - predicted_deaths,
    residual_deaths_pe =
      (observed_deaths - predicted_deaths) / observed_deaths * 100,
    residual_mortality_pe =
      (observed_mortality - predicted_mortality) / observed_mortality * 100,
    residual_log_mortality_e =
      log(observed_mortality) - log(predicted_mortality)
  )

# Residual summaries by model and sex -----------------------------

dat$summarised_residuals <-
  dat$residuals_by_week_and_sex %>%
  group_by(model_name, sample, sex) %>%
  summarise(
    mpe = mean(residual_deaths_pe, na.rm = TRUE),
    me = mean(residual_deaths_e, na.rm = TRUE),
    mdape = median(abs(residual_deaths_pe), na.rm = TRUE),
    mdae = median(abs(residual_deaths_e), na.rm = TRUE)
  ) %>%
  ungroup()

# Plot test errors over week --------------------------------------

dat$residuals_by_week_and_sex %>%
  pivot_longer(cols = c(residual_deaths_e:residual_log_mortality_e),
               names_to = 'residual_type', values_to = 'residual_value') %>%
  filter(sample == 'test') %>%
  group_by(residual_type) %>%
  group_walk(~{
    fig[[paste0(.y$residual_type, '_by_week')]] <<-
      .x %>%
      ggplot(aes(x = iso_week, y = residual_value)) +
      geom_point(
        color = 'grey60', size = 0.3
      ) +
      geom_smooth(
        se = FALSE, color = 'black',
      ) +
      geom_hline(yintercept = 0) +
      glob$ggtheme +
      facet_grid(sex~model_name) +
      scale_x_continuous(breaks = c(1, seq(10, 40, 10), 52)) +
      labs(x = 'Week of year', y = paste0('Residual ', .y$residual_type))
  })

# distribution of residuals on test set
dat$residuals_by_week_and_sex %>%
  pivot_longer(cols = c(residual_deaths_e:residual_log_mortality_e),
               names_to = 'residual_type', values_to = 'residual_value') %>%
  filter(sample == 'test') %>%
  group_by(residual_type) %>%
  group_walk(~{
    residual_summary <-
      .x %>%
      group_by(model_name, sex) %>%
      summarise(
        mean = mean(residual_value, na.rm = TRUE),
        median = median(residual_value, na.rm = TRUE)
      ) %>%
      ungroup()
    
    fig[[paste0(.y$residual_type, '_distribution')]] <<-
      ggplot(.x) +
      geom_vline(xintercept = 0, color = 'grey80') +
      geom_density(
        aes(
          x = residual_value,
          color = model_name
        )
      ) +
      geom_text(
        aes(
          x = 0, y = 0,
          label = paste0('Mean ', formatC(mean, 3)),
          color = model_name
        ),
        size = 3,
        data = residual_summary,
        hjust = 0.5,
        position = position_stack(vjust = 1)
      ) +
      scale_x_continuous() +
      facet_grid(sex~model_name) +
      guides(color = 'none') +
      glob$MyGGplotTheme() +
      labs(x = .y$residual_type)
  })

# Plot test accuracy by model -------------------------------------

dat$summarised_residuals %>%
  pivot_longer(
    c(mpe, me, mdape, mdae),
    names_to = 'summary_type', values_to = 'summary_value'
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = summary_value
  ) %>%
  group_by(summary_type) %>%
  group_walk(~{
    fig[[paste0(.y$summary_type, '_summary')]] <<-
      .x %>%
      mutate(
        model_name =
          fct_reorder(model_name, test, .desc = FALSE, .fun = mean)
      ) %>%
      ggplot(aes(group = sex, color = sex)) +
      geom_linerange(
        aes(x = model_name, ymin = training, ymax = test),
        #arrow = arrow(angle = 10, type = 'closed', length = unit(3, 'mm')),
        position = position_dodge2(width = 0.5)
      ) +
      geom_point(
        aes(x = model_name, y = training),
        position = position_dodge2(width = 0.5),
        fill = 'white',
        shape = 21
      ) +
      geom_point(
        aes(x = model_name, y = test),
        position = position_dodge2(width = 0.5)
      ) +
      geom_text(
        aes(x = model_name, y = test+0.4,
            label = paste0(formatC(test, digits = 3))),
        position = position_dodge2(width = 0.5)
      ) +
      glob$MyGGplotTheme() +
      coord_flip() +
      labs(x = NULL, y = .y$summary_type)
    
  })

# Export plots ----------------------------------------------------

ExportFiguresFromList(fig, path = 'out', add_date = TRUE)
