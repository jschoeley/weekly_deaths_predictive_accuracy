# Weekly mortality/deaths predictive accuracy
# Jonas Schöley

# Init ------------------------------------------------------------

library(tidyverse)
library(mgcv)

cnst <- list(
  country = 'DNK',
  week_season_year_starts = 27,
  ggtheme =
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(linetype = 3, color = "grey80")
    ),
  col_training =
    c(training = 'grey70', test = 'grey30')
)

dat <- list()
fig <- list()
tab <- list()

# Functions -------------------------------------------------------

#' Convert Week of Year to Date
#'
#' @param year Year integer.
#' @param week Week of year integer (1 to 53).
#' @param weekday Weekday integer (1, Monday to 7, Sunday).
#' @param offset Integer offset added to `week` before date calculation.
#'
#' @return A date object.
#' 
#' @source https://en.wikipedia.org/wiki/ISO_8601
#'
#' @author Jonas Schöley
#'
#' @examples
#' # the first Week of 2020 actually starts Monday, December 30th 2019
#' ISOWeekDate2Date(2020, 1, 1)
ISOWeekDate2Date <- function (year, week, weekday = 1, offset = 0) {
  require(ISOweek)
  isoweek_string <-
    paste0(
      year, '-W',
      formatC(
        week+offset,
        flag = '0',
        format = 'd',
        digits = 1
      ),
      '-', weekday
    )
  ISOweek2date(isoweek_string)
}

#' Calculate Weeks Since Some Origin Date
#'
#' @param date Date string.
#' @param origin_date Date string.
#' @param week_format Either 'integer' for completed weeks or
#' 'fractional' for completed fractional weeks.
#'
#' @return Time difference in weeks.
#' 
#' @author Jonas Schöley
#'
#' @examples
#' # My age in completed weeks
#' WeeksSinceOrigin('2020-07-07', '1987-07-03')
WeeksSinceOrigin <-
  function (date, origin_date, week_format = 'integer') {
    require(ISOweek)
    fractional_weeks_since_origin <-
      as.double(difftime(
        as.Date(date),
        as.Date(origin_date),
        units = 'weeks'
      ))
    switch(
      week_format,
      fractional = fractional_weeks_since_origin,
      integer = as.integer(fractional_weeks_since_origin)
    )
  }

SeasonYearSequence <- function (from, to) {
  years <- from:to
  paste0(head(years, -1), '/', years[-1])
}

# This model and estimates the average mortality rate over
# some years within each week and stratum. The associated
# predict() method multiplies this average mortality with
# given exposures to derive death counts.

AverageMortalityModel <-
  function (df, week, deaths, exposures, ...) {
    require(dplyr)
    .strata = enquos(...); .week = enquo(week);
    .deaths = enquo(deaths); .exposures = enquo(exposures)
    
    avg_mx <-
      df %>%
      group_by(!!!.strata, !!.week) %>%
      summarise(
        avg_mortality = mean(!!.deaths/!!.exposures),
        .groups = 'drop'
      )
    
    structure(list(avg = avg_mx), class = 'avgmx')
    
  }

predict.avgmx <- function (object, newdata, ...) {
  require(dplyr)
  
  suppressMessages(left_join(newdata, object$avg)) %>%
    pull(avg_mortality)
}

ExportPNG <-
  function(figure,
           filename,
           path,
           width = 170,
           height = 100,
           scale = 1,
           ...) {
    require(ggplot2)
    ggsave(
      filename = paste0(filename, ".png"),
      plot = figure,
      path = path,
      width = width,
      height = height,
      units = "mm",
      dpi = 300,
      scale = scale,
      ...
    )
  }

ExportPDF <-
  function(figure,
           filename,
           path,
           width = 170,
           height = 100,
           scale = 1,
           ...) {
    require(ggplot2)
    ggsave(
      filename = paste0(filename, ".pdf"),
      plot = figure,
      path = path,
      width = width,
      height = height,
      units = "mm",
      dpi = 300,
      scale = scale,
      useDingbats = FALSE,
      ...
    )
  }

# Load data -------------------------------------------------------

dat$raw <-
  readr::read_csv(
    "data/stmf-2.csv",
    col_types = "ciicddddddddddddlll",
    skip = 2, col_names = TRUE
  )

# Prepare data ----------------------------------------------------

dat$ts <-
  dat$raw %>%
  filter(
    CountryCode == cnst$country,
    Sex != 'b'
  ) %>%
  select(
    year = Year, week = Week, sex = Sex, D0_14:D85p, R0_14:R85p,
  ) %>%
  pivot_longer(
    cols = D0_14:R85p,
    names_sep = 1,
    names_to = c('statistic', 'age_group'),
    values_to = 'value'
  ) %>%
  pivot_wider(
    names_from = statistic,
    values_from = value
  ) %>%
  # discard lowest age group as it carries almost no information,
  # yet clutters the display of results
  #filter(age_group != '0_14') %>%
  rename(observed_deaths = D, observed_mortality = R) %>%
  mutate(
    # use same weekly exposure as used by HMD
    exposure =
      observed_deaths / observed_mortality,
    date =
      ISOWeekDate2Date(year, week, 1),
    start_of_season_year =
      ISOWeekDate2Date(
        ifelse(week<=cnst$week_season_year_starts, year-1, year),
        cnst$week_season_year_starts, 1
      ),
    season_year =
      paste0(lubridate::year(start_of_season_year),'/',
             lubridate::year(start_of_season_year)+1),
    season_year_int =
      as.integer(substr(season_year,1,4)),
    weeks_into_season_year =
      difftime(date, start_of_season_year, units = 'weeks') %>%
      floor() %>% as.integer()
  )

# Create training/test sets ---------------------------------------

dat$ts %>%
  mutate(m = observed_deaths/exposure) %>%
  ggplot() +
  geom_point(aes(x = date, y = m)) +
  scale_x_date(date_breaks = '3 months', date_labels = '%b') +
  facet_wrap(~season_year, scales = 'free_x')

dat$tt <-
  # set up K-fold cross-validation
  #SWE map2(2000+0:12, 2004+0:12, SeasonYearSequence) %>%
  #ENW map2(2010+0:5, 2014+0:5, SeasonYearSequence) %>%
  map2(2006+0:8, 2011+0:8, SeasonYearSequence) %>%
  map(~filter(dat$ts, season_year %in% .x)) %>%
  map(~mutate(
    .x, training = ifelse(
      season_year_int == max(season_year_int) |
        (year == max(season_year_int) & week >= 10),
      'test', 'training')
  )) %>%
  bind_rows(.id = 'cv_id') %>%
  mutate(cv_id = as.integer(cv_id))

# Plot training-test split ------------------------------------------------

# plot training-test split
fig$training_test_bars <-
  dat$tt %>%
  filter(age_group == '85p', sex == 'm') %>%
  ggplot(aes(x = date, y = cv_id)) +
  geom_path(aes(color = training, group = interaction(cv_id, training)), size = 3) +
  cnst$ggtheme +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  scale_y_continuous(breaks = 1:20) +
  labs(x = NULL, y = NULL) +
  guides(color = 'none') +
  scale_color_manual(values = cnst$col_training)

fig$weekly_observed_deaths <-
  dat$tt %>%
  group_by(date, cv_id) %>%
  summarise(
    observed_deaths = sum(observed_deaths)
  ) %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = observed_deaths)) +
  cnst$ggtheme +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  labs(x = NULL, y = 'Weekly observed deaths')

ExportPDF(fig$weekly_observed_deaths,
          'weekly_observed_deaths',
          path = './out')

fig$cross_validation <-
  dat$tt %>%
  group_by(date, cv_id, training) %>%
  summarise(
    observed_deaths = sum(observed_deaths)
  ) %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = cv_id + 0.001*observed_deaths-1,
                 color = training), size = 0.1) +
  geom_point(aes(y = 0.001*observed_deaths-1), size = 0.1) +
  cnst$ggtheme +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  scale_y_continuous(breaks = 1:20) +
  labs(x = NULL, y = 'Cross Validation Series') +
  guides(color = 'none') +
  scale_color_manual(values = cnst$col_training)

ExportPDF(fig$cross_validation,
          'cross_validation',
          path = './out')

# Specifications of models to test ----------------------------------------

dat$mod_spec <-
  tribble(
    ~model_name, ~model_class, ~model_spec,
    'Avg. weekly mortality', 'avgmx', NA,
    'Serfling-Poisson GLM', 'glm', list(
      formula = formula(
        observed_deaths ~
          # log linear long term trend
          weeks_since_origin*sex*age_group +
          # seasonality
          # full year period
          sin(2*pi*weeks_into_season_year/(365.25/7))*sex*age_group +
          cos(2*pi*weeks_into_season_year/(365.25/7))*sex*age_group +
          # half year period
          sin(2*pi*weeks_into_season_year/(365.25/2/7))*sex*age_group +
          cos(2*pi*weeks_into_season_year/(365.25/2/7))*sex*age_group +
          # adjustment for new years eve
          new_year*sex*age_group +
          # exposures
          offset(log(exposure))
      ),
      family = quasipoisson(link = 'log')
    ),
    'Poisson GAM', 'gam', list(
      formula = formula(
        observed_deaths ~
          1 + sex + age_group +
          # log linear long term trend
          weeks_since_origin*sex_age_interaction +
          # penalized cyclic spline for seasonality
          s(weeks_into_season_year, bs = 'cp', k = 52, by = sex_age_interaction) +
          # adjustment for new years eve
          new_year*sex_age_interaction +
          # exposures
          offset(log(exposure))
      ),
      family = quasipoisson(link = 'log')
    ),
    'Neg-Bin. GAM', 'gam', list(
      formula = formula(
        observed_deaths ~
          1 + sex + age_group +
          # log linear long term trend
          weeks_since_origin*sex_age_interaction +
          # penalized cyclic spline for seasonality
          s(weeks_into_season_year, bs = 'cp', k = 52, by = sex_age_interaction) +
          # adjustment for new years eve
          new_year*sex_age_interaction +
          # exposures
          offset(log(exposure))
      ),
      family = nb(link = 'log')
    ),
    'Neg-Bin GAM-RE', 'gam', list(
      formula = formula(
        observed_deaths ~
          1 + sex + age_group +
          # log linear long term trend
          weeks_since_origin*sex_age_interaction +
          # penalized cyclic spline for seasonality
          s(weeks_into_season_year, bs = 'cp', k = 52, by = sex_age_interaction) +
          # adjustment for new years eve
          new_year*sex_age_interaction +
          # season year random effect adjustment
          s(season_year_sex_age_interaction, bs = 're') +
          # exposures
          offset(log(exposure))
      ),
      family = nb(link = 'log')
    )
  ) %>%
  mutate(model_id = 1:n())

# Model time series -----------------------------------------------

dat$tt <-
  dat$tt %>%
  nest(
    data = c(-cv_id)
  ) %>%
  expand_grid(
    nest(dat$mod_spec, model_spec = c(-model_name, -model_id, -model_class))
  )

dat$mod <-
  dat$tt %>%
  group_by(model_name, cv_id) %>%
  group_modify(~{
    
    cat(unlist(.y), sep = '\n')
    
    ### PREPARE INPUT DATA
    
    # season-years in training data
    # I use this later to check if prediction
    # outside the range of season_years in training
    # data will be required
    season_years_in_training <-
      .x$data[[1]] %>%
      filter(training == 'training') %>%
      pull(season_year) %>% unique()
    
    # prepare for model
    model_data <-
      .x$data[[1]] %>%
      mutate(
        weeks_since_origin =
          WeeksSinceOrigin(date, min(date)),
        observed_mortality =
          observed_deaths / exposure,
        new_year =
          ifelse(week == 52, 1, 0),
        sex_age_interaction =
          interaction(sex, age_group),
        include_re =
          ifelse(season_year %in% season_years_in_training, 1, 0),
        season_year_sex_age_interaction =
          interaction(season_year, sex, age_group)
      )
    
    training_data <-
      model_data %>%
      filter(training == 'training') %>%
      mutate(weight = exp(weeks_since_origin/max(weeks_since_origin)))
    
    ### TRAIN MODELS
    
    if (.x$model_class == 'avgmx') {
      
      # calculate average weekly mortality rate over training data
      model <-
        AverageMortalityModel(
          training_data,
          week,
          observed_deaths,
          exposure,
          sex, age_group
        )
      
    }
    
    if (.x$model_class == 'glm') {
      
      model <-
        glm(
          formula = .x$model_spec[[1]][[1]][[1]]$formula,
          data = training_data,
          family = .x$model_spec[[1]][[1]][[1]]$family
        )
      
    }
    
    if (.x$model_class == 'gam') {
      
      model <-
        gam(
          formula = .x$model_spec[[1]][[1]][[1]]$formula,
          data = training_data,
          family = .x$model_spec[[1]][[1]][[1]]$family,
          method = 'REML'
        )

    }

    ### PREDICT FROM MODEL
    
    predictions <-
      model_data %>%
      group_by(include_re) %>%
      group_modify(~{
        newdata <- mutate(.x, exposure = 1)
        # this is to allow predictions with novel levels in
        # random effects variables. for the novel levels the
        # random effect is then set to 0.
        exclude = NULL
        if (.y == 0) {
          newdata$season_year_sex_age_interaction <- NULL
          exclude <- 's(season_year_sex_age_interaction)'
        }
        # predict mortality rates
        mutate(
          .x,
          predicted_mortality =
            predict(
              model,
              newdata = newdata,
              type = 'response',
              # for prediction with random effects
              exclude = exclude,
              newdata.guaranteed = TRUE
            ),
          predicted_deaths =
            predicted_mortality * exposure)
      })
    
    result <-
      tibble(
        training_data = list(training_data),
        predictions = list(predictions),
        fitted_model = list(model)
      )
    
    return(result)
    
  }) %>%
  ungroup()

# mod <- dat$mod
#save(mod, file = './out/model_fits.RData')
#load('./out/model_fits.RData')
#dat$mod <- mod

# Derive residuals --------------------------------------------------------

dat$residuals_by_week_and_sex <-
  dat$mod %>%
  unnest(predictions) %>%
  group_by(cv_id, model_name, training,
           week, sex) %>%
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
      (observed_mortality - predicted_mortality) / observed_mortality * 100
  )

dat$summarised_residuals_by_week_and_sex <-
  dat$residuals_by_week_and_sex %>%
  group_by(model_name, training, sex) %>%
  summarise(
    mdape = median(abs(residual_deaths_pe), na.rm = TRUE),
    mdae = median(abs(residual_deaths_e), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from = training,
    values_from = c(mdape, mdae),
    id_cols = c(model_name, sex)
  ) %>%
  mutate(
    model_name =
      fct_reorder(model_name, mdape_test, .desc = TRUE, .fun = mean)
  )

# Plot observed vs. fitted ------------------------------------------------

dat$mod_spec %>%
  group_by(model_name) %>%
  group_walk(~{
    scaler <- 0.002
    fig_dat <-
      dat$mod %>%
      filter(model_name == .y[[1]]) %>%
      unnest(predictions) %>%
      group_by(model_name, cv_id, date, training) %>%
      summarise(
        observed_mortality = sum(observed_deaths, na.rm = T),
        predicted_mortality = sum(predicted_deaths, na.rm = T)
      )
    
    dat$fig[[.y[[1]]]] <<-
      ggplot(fig_dat, aes(x = date)) +
      geom_point(
        aes(
          group = cv_id,
          color = training,
          y = cv_id + scaler*observed_mortality-2
        ),
        size = 0.3
      ) +
      geom_point(aes(y = scaler*observed_mortality-2), size = 0.3) +
      geom_line(aes(y = cv_id + scaler*predicted_mortality-2,
                    group = cv_id, alpha = training),
                color = 'red') +
      facet_wrap(~model_name, nrow = 2) +
      cnst$ggtheme +
      scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
      scale_y_continuous(breaks = 1:9) +
      scale_alpha_manual(values = c(training = 0.3, test = 1)) +
      labs(x = NULL, y = 'Cross Validation Series') +
      guides(color = 'none', alpha = 'none') +
      scale_color_manual(values = cnst$col_training)
  })

ExportPNG(dat$fig$`Avg. weekly mortality`,
          'avg_weekly_mortality',
          path = './out')

ExportPNG(dat$fig$`Serfling-Poisson GLM`,
          'serfling_poisson_glm',
          path = './out')

ExportPNG(dat$fig$`Poisson GAM`,
          'poisson_gam',
          path = './out')

ExportPNG(dat$fig$`Neg-Bin. GAM`,
          'negbin_gam',
          path = './out')

ExportPNG(dat$fig$`Neg-Bin GAM-RE`,
          'negbin_gam_re',
          path = './out')

# Plot residuals ----------------------------------------------------------

dat$fig$residual_deaths_pe <-
  dat$residuals_by_week_and_sex %>%
  ggplot(aes(x = week, y = residual_deaths_pe)) +
  geom_point(
    aes(color = training), alpha = 0.5
  ) +
  geom_smooth(
    aes(group = training, lty = training),
    se = FALSE, color = 'black',
  ) +
  geom_hline(yintercept = 0) +
  cnst$ggtheme +
  facet_grid(sex~model_name) +
  scale_x_continuous(breaks = c(1, seq(10, 40, 10), 52)) +
  labs(x = 'Week of year', y = 'Residual deaths percent error')

ExportPDF(dat$fig$residual_deaths_pe, 'residual_deaths_pe', './out',
          scale = 1.3)

dat$fig$residual_deaths_e <-
  dat$residuals_by_week_and_sex %>%
  ggplot(aes(x = week, y = residual_deaths_e)) +
  geom_point(
    aes(color = training), alpha = 0.5
  ) +
  geom_smooth(
    aes(group = training, lty = training),
    se = FALSE, color = 'black',
  ) +
  geom_hline(yintercept = 0) +
  cnst$ggtheme +
  facet_grid(sex~model_name) +
  scale_x_continuous(breaks = c(1, seq(10, 40, 10), 52)) +
  labs(x = 'Week of year', y = 'Residual deaths')

ExportPDF(dat$fig$residual_deaths_e, 'residual_deaths_e', './out',
          scale = 1.3)

dat$residuals_by_week_and_sex %>%
  filter(training == 'test') %>%
  ggplot() +
  geom_density(
    aes(
      x = residual_deaths_e,
      color = model_name
    )
  ) +
  scale_x_continuous() +
  geom_vline(xintercept = 0) +
  facet_wrap(~sex)

# Plot predictive accuracy ------------------------------------------------

dat$fig$mdape <-
  dat$summarised_residuals_by_week_and_sex %>%
  ggplot(aes(group = sex, color = sex)) +
  geom_linerange(
    aes(x = model_name, ymin = mdape_training, ymax = mdape_test),
    #arrow = arrow(angle = 10, type = 'closed', length = unit(3, 'mm')),
    position = position_dodge2(width = 0.5)
  ) +
  geom_point(
    aes(x = model_name, y = mdape_training),
    position = position_dodge2(width = 0.5),
    fill = 'white',
    shape = 21
  ) +
  geom_point(
    aes(x = model_name, y = mdape_test),
    position = position_dodge2(width = 0.5)
  ) +
  geom_text(
    aes(x = model_name, y = mdape_test+0.4,
        label = paste0(formatC(mdape_test, digits = 2, format = 'f'))),
    position = position_dodge2(width = 0.5)
  ) +
  cnst$ggtheme +
  coord_flip() +
  labs(x = NULL, y = 'Median absolute percentage error')

ExportPNG(dat$fig$mdape,
          'mdape',
          path = './out')

# female weekly death counts are harder to predict than males
# average mortality method is positively biased as mortality improvements matter
# also in the short term
# little difference between different log-linear regression model specifications
