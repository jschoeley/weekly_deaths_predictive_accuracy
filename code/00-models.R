# Short-term mortality forecasting models
# 
# 2020-10-19
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
  formula,
  winter_deaths_epi_weeks, exogenous_epi_weeks, reference_week,
  stratum_name, week_name, death_name, exposure_name, sample_name, year_name
) {
  
  require(dplyr)
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
    ungroup()
  
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
      # important for the alr transformation
      deaths_reference_week_yj =
        sum(ifelse(!!.week %in% reference_week, !!.deaths, 0)),
      share_on_annual_deaths_reference_week_yj =
        deaths_reference_week_yj / annual_deaths_yj,
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
    # alr transform
    dat_alr <-
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
        # outcome variable
        share_on_annual_deaths_wyj =
          !!.deaths / annual_deaths_yj,
        alr_share_on_annual_deaths =
          log(share_on_annual_deaths_wyj+fudge) -
          log(share_on_annual_deaths_reference_week_yj+fudge),
        alr_share_on_annual_deaths =
          ifelse(!!.week == reference_week, NA, alr_share_on_annual_deaths),
        # predictors
        deaths_winter_centered_yj =
          (deaths_winter_yj - avg_deaths_winter_j) / sd_deaths_winter_j,
        mortality_winter_centered_yj =
          (mortality_winter_yj - avg_mortality_winter_j) / sd_mortality_winter_j
      ) %>%
      ungroup()
    
    # remove reference week
    dat_alr_no_reference_week <-
      dat_alr %>%
      filter(!!.week != reference_week) %>%
      mutate(w_fac = factor(!!.week))
    
    # prepare training and prediction data
    dat_train <- filter(dat_alr_no_reference_week, epi_year_contains_test_data == FALSE)
    dat_pred <- dat_alr_no_reference_week
    
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
        # predicted alr transformed death proportions by
        # epi-year, sex, age and week
        pred_alr_yjw =
          predict(model, newdata = ., type = 'response')
      ) %>%
      # merge into original data which does include reference week
      right_join(dat_alr) %>%
      # transform alr back to proportions and then deaths
      mutate(
        pred_exp_alr_yjw =
          ifelse(!!.week == reference_week, 1, exp(pred_alr_yjw))
      ) %>%
      group_by(!!.year, !!.stratum_id) %>%
      mutate(
        pred_share_on_annual_deaths_yjw =
          pred_exp_alr_yjw/sum(pred_exp_alr_yjw),
        pred_share_exogenous_deaths_yj =
          sum(ifelse(!!.week %in% exogenous_epi_weeks, pred_share_on_annual_deaths_yjw, 0))
      ) %>%
      ungroup() %>%
      mutate(
        pred_annual_deaths_yj =
          deaths_exogenous_yj/pred_share_exogenous_deaths_yj,
        predicted_deaths =
          pred_share_on_annual_deaths_yjw*pred_annual_deaths_yj
      )
  })

  #predictions <- df_input
  #predictions$predicted_deaths <- dat_pred$predicted_deaths
  
  return(
    list(model = model, predictions = predictions)
  )
  
}
