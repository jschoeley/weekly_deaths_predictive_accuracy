# Short-term mortality forecasting models
# 
# 2020-10-20
#
# Jonas Schöley

# All models take a data frame with training data
# which is used for model estimation, called <df_training>,
# and a data frame used to make prediction given the fitted model,
# called <df_prediction>.
#
# All models output a list with two elements called <model>,
# an object containing the fitted model, coefficients etc., and
# <predictions>, which is the input data frame <df_prediction> with
# an added column called <predicted_deaths>.

# For test purposes -----------------------------------------------

load('data/2020-10-05-xstmf_cv.RData')
library(dplyr)
test_data <-
  xstmf_cv %>%
  filter(country_code == 'DNK', cv_id == 1) %>%
  rename(y = epi_year, w = epi_week, j = stratum_id) %>%
  mutate(
    j_fac = factor(j),
    w_fac = as.factor(w),
    y_short = paste0(substr(y, 3, 4), '/', substr(y, 8, 9))
  )

df_input <- test_data
df_training <- filter(test_data, sample == 'training')
df_test <- filter(test_data, sample == 'test')

.week <- quo(w)
.deaths <- quo(observed_deaths)
.exposures <- quo(exposure_pw_hmd)
.year <- quo(y)
.sample <- quo(sample)
.stratum_id <- quo(j_fac)

# Simple averages model -------------------------------------------

# This model and estimates the average mortality rate over
# some years within each week and stratum. The associated
# predict() method multiplies this average mortality with
# given exposures to derive death counts.
SimpleAveragesModel <- function (
  df_training, df_prediction,
  week_name, death_name, exposure_name, year_name,
  n_years,
  stat = 'mortality', ...
) {
  
  require(dplyr)
  .strata <- enquos(...)
  .week <- enquo(week_name)
  .deaths <- enquo(death_name)
  .exposures <- enquo(exposure_name)
  .year <- enquo(year_name)

  # last <n_years> number of years before start of test data
  training_years <-
    df_training %>% pull(!!.year) %>% unique() %>%
    sort(decreasing = TRUE) %>% `[`(1:n_years)
    
  df_training_grouped_filtered <-
    df_training %>%
    filter(!!.year %in% training_years) %>%
    group_by(!!!.strata, !!.week)
    
  if (identical(stat, 'mortality')) {
    average <-
      df_training_grouped_filtered %>%
      summarise(avg = mean(!!.deaths / !!.exposures), .groups = 'drop')
    predictions <-
      suppressMessages(left_join(df_prediction, average)) %>%
      mutate(predicted_deaths = avg*!!.exposures) %>%
      select(-avg)
  }
  
  if (identical(stat, 'deaths')) {
    average <-
      df_training_grouped_filtered %>%
      summarise(avg = mean(!!.deaths), .groups = 'drop')
    predictions <-
      suppressMessages(left_join(df_prediction, average)) %>%
      mutate(predicted_deaths = avg) %>%
      select(-avg)
  }
  
  model <- structure(
    list(
      averages = average,
      strata_name = .strata,
      week_name = .week,
      death_name = .deaths,
      exposure_name = .exposures
    ),
    class = 'avg'
  )
  
  return(
    list(model = model, predictions = predictions)
  )
  
}

# Serfling model --------------------------------------------------

SerflingModel <- function (
  df_training, df_prediction,
  formula, family
) {
  
  require(dplyr)
  
  model <- glm(
    formula = formula,
    family = family,
    data = df_training
  )
  
  predictions <-
    df_prediction %>%
    mutate(
      predicted_deaths =
        predict(model, newdata = ., type = 'response')
    )
  
  return(
    list(model = model, predictions = predictions)
  )
  
}

# GAM model -------------------------------------------------------

GAMModel <- function (
  df_training, df_prediction,
  formula, family, method
) {
  
  require(dplyr)
  require(mgcv)
  
  model <- gam(
    formula = formula,
    family = family,
    data = df_training,
    method = method
  )
  
  predictions <-
    df_prediction %>%
    mutate(
      predicted_deaths =
        predict(model, newdata = ., type = 'response')
    )
  
  return(
    list(model = model, predictions = predictions)
  )
  
}

# CODA model ------------------------------------------------------

CODAModel <- function (
  df_input,
  formula, transform,
  winter_deaths_epi_weeks, exogenous_epi_weeks,
  stratum_name, week_name, death_name, exposure_name, sample_name, year_name
) {
  
  require(dplyr)
  require(compositions)
  .stratum_id <- enquo(stratum_name)
  .week <- enquo(week_name)
  .deaths <- enquo(death_name)
  .exposures <- enquo(exposure_name)
  .sample <- enquo(sample_name)
  .year <- enquo(year_name)
  
  # constants
  fudge = 1e-6 # yeah, no 0s allowed in CODA...
  
  # because coda can't train on partial information for a year
  # add an indicator stating if a given year contains test data
  dat_input <-
    df_input %>%
    group_by(!!.stratum_id, !!.year) %>%
    mutate(epi_year_contains_test_data = any(!!.sample == 'test')) %>%
    ungroup() %>%
    mutate(
      # add row id for later identification
      row_id = 1:n()
    )
  
  # year and stratum specific summary variables
  dat_yj <-
    dat_input %>%
    group_by(!!.stratum_id, !!.year) %>%
    summarise(
      # predictors
      deaths_winter_yj =
        sum(ifelse(!!.week %in% winter_deaths_epi_weeks, !!.deaths, 0)),
      exposure_winter_yj =
        sum(ifelse(!!.week %in% winter_deaths_epi_weeks, !!.exposures, 0)),
      mortality_winter_yj =
        deaths_winter_yj / exposure_winter_yj,
      # for deriving total deaths later on
      deaths_exogenous_yj =
        sum(ifelse(!!.week %in% exogenous_epi_weeks, !!.deaths, 0)),
      annual_deaths_yj =
        # annual deaths are unknown for years where only partial
        # training data is available
        sum(ifelse(epi_year_contains_test_data, NA, !!.deaths)),
      .groups = 'drop'
    )
  
  # stratum specific summary variables
  dat_j <-
    dat_yj %>%
    group_by(!!.stratum_id) %>%
    summarise(
      avg_deaths_winter_j =
        mean(deaths_winter_yj),
      sd_deaths_winter_j =
        sd(deaths_winter_yj),
      avg_mortality_winter_j =
        mean(mortality_winter_yj),
      sd_mortality_winter_j =
        sd(mortality_winter_yj),
      .groups = 'drop'
    )
  suppressMessages({
    # add multilevel predictors to input data and compute
    # share on annual death by week per stratum
    dat_jyw <-
      dat_input %>%
      left_join(
        dat_yj,
        by = c(as_label(.stratum_id), as_label(.year))
      ) %>%
      left_join(
        dat_j,
        by = c(as_label(.stratum_id))
      ) %>%
      mutate(
        # predictors
        deaths_winter_centered_yj =
          (deaths_winter_yj - avg_deaths_winter_j) / sd_deaths_winter_j,
        mortality_winter_centered_yj =
          (mortality_winter_yj - avg_mortality_winter_j) / sd_mortality_winter_j,
        # outcome variable
        share_on_annual_deaths_wyj =
          !!.deaths / annual_deaths_yj
      ) %>%
      ungroup()
    
    # transform outcome variable to matrix with one column
    # per week and one row per case (i.e., stratum and year)
    dat_jyw <-
      dat_jyw %>%
      arrange(!!.stratum_id, !!.year, !!.week) %>%
      group_by(!!.stratum_id, !!.year) %>%
      group_modify(.f = ~{
        x <- .x[['share_on_annual_deaths_wyj']]
        # number of weeks in composition
        # may be 52 or 53
        n_weeks <- length(pull(.x, !!.week))
        if (any(is.na(x))) {
          # if part of the annual share is missing
          # don't bother with a CODA transform
          z <- NA
        } else {
          z <- switch(
            transform,
            ilr = ilr(x+fudge, V = ilrBase(D=n_weeks)),
            alr = alr(x+fudge, ivar = n_weeks) 
          )
          z <- c(z, NA)
        }
        mutate(
          .x, transformed_share_on_annual_deaths = z,
          n_weeks = n_weeks
        )
      }) %>%
      ungroup()
    
    dat_jyw_no_reference_week <-
      dat_jyw %>%
      # within each stratum-year remove the reference week
      group_by(!!.stratum_id, !!.year) %>%
      slice(-n_weeks[1]) %>%
      mutate(across(!!.week, factor)) %>%
      ungroup()
    # prepare training and prediction data
    dat_train <- filter(dat_jyw_no_reference_week, epi_year_contains_test_data == FALSE)
    dat_pred <- dat_jyw_no_reference_week
    
    # fit model
    model <-
      glm(
        formula = formula,
        family = gaussian(link = 'identity'),
        data = dat_train
      )
    
    # predict from fitted model
    predictions <-
      dat_pred %>%
      mutate(
        # predicted transformed weekly death proportions by
        # epi-year and stratum
        y_hat =
          predict(model, newdata = ., type = 'response')
      ) %>%
      select(row_id, y_hat) %>%
      # merge into original data which does include reference week
      right_join(dat_jyw, by = 'row_id') %>%
      # within each stratum-year, convert predicted transformed
      # weekly death proportions to proportion scale
      group_by(!!.stratum_id, !!.year) %>%
      group_modify(~{
        # predicted transformed share without reference week
        z_hat <- .x$y_hat[-(.x$n_weeks[1])]
        # back-transform from coda space to proportion space
        x_hat <- switch(
          transform,
          ilr = ilrInv(z_hat, V = ilrBase(D=.x$n_weeks[1])),
          alr = alrInv(z_hat)
        )
      mutate(.x, pred_share_on_annual_deaths_yjw = x_hat)
      }) %>%
      mutate(
        pred_share_exogenous_deaths_yj =
          sum(ifelse(!!.week %in% exogenous_epi_weeks, pred_share_on_annual_deaths_yjw, 0))
      ) %>%
      select(-n_weeks) %>%
      ungroup() %>%
      # convert to death count scale
      mutate(
        pred_annual_deaths_yj =
          deaths_exogenous_yj/pred_share_exogenous_deaths_yj,
        predicted_deaths =
          pred_share_on_annual_deaths_yjw*pred_annual_deaths_yj
      ) %>%
      select(row_id, predicted_deaths) %>%
      right_join(dat_input, by = 'row_id') %>%
      select(-row_id)
  })
  
  return(
    list(model = model, predictions = predictions)
  )
  
}
