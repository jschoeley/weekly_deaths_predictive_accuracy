# Predict weekly death counts by sex and age
#
# 2020-10-15
#
# Jonas Schöley

# Init ------------------------------------------------------------

set.seed(1987)
source('code/00-init.R')
source('code/00-models.R')

library(tidyverse)

cnst <- list(
  # which countries to analyze?
  countries = c('DNK', 'SWE')
)

dat <- list()
fig <- list()

# Load data -------------------------------------------------------

load('data/2020-10-05-xstmf_cv.RData')
dat$cv_data <-
  xstmf_cv %>%
  filter(country_code %in% cnst$countries) %>%
  rename(y = epi_year, w = epi_week, j = stratum_id) %>%
  mutate(
    j_fac = factor(j),
    w_fac = as.factor(w),
    y_short = paste0(substr(y, 3, 4), '/', substr(y, 8, 9))
  )

# Specifications of models to test --------------------------------

# Just a big list specifying all the models to be tested,
# and their parametrizations. See 00-models.R for the exact
# implementation of the models.

dat$mod_spec <-
  tribble(
    ~model_name, ~model_name_short, ~model_class, ~model_spec,
    '5-year avg. weekly mortality', '5y-avg-mort', 'avg', list(
      stat = 'mortality',
      n_years = 5
    ),
    '3-year avg. weekly mortality', '3y-avg-mort', 'avg', list(
      stat = 'mortality',
      n_years = 3
    ),
    '5-year avg. weekly deaths', '5y-avg-death', 'avg', list(
      stat = 'deaths',
      n_years = 5
    ),
    '3-year avg. weekly deaths', '3y-avg-death', 'avg', list(
      stat = 'deaths',
      n_years = 3
    ),
    'Serfling-Poisson GLM', 'serfl', 'serfling', list(
      formula = formula(
        observed_deaths ~
          # log linear long term trend
          weeks_since_origin*j +
          # seasonality
          # full year period
          sin(2*pi*w/(365.25/7))*j_fac +
          cos(2*pi*w/(365.25/7))*j_fac +
          # half year period
          sin(2*pi*w/(365.25/2/7))*j_fac +
          cos(2*pi*w/(365.25/2/7))*j_fac +
          # adjustment for special weeks
          special_week*j_fac +
          # exposures
          offset(log(exposure_pw_hmd))
      ),
      family = quasipoisson(link = 'log')
    ),
    'Serfling-Poisson GLM (no exposures)', 'serfl-no_expos', 'serfling', list(
      formula = formula(
        observed_deaths ~
          # log linear long term trend
          weeks_since_origin*j_fac +
          # seasonality
          # full year period
          sin(2*pi*w/(365.25/7))*j_fac +
          cos(2*pi*w/(365.25/7))*j_fac +
          # half year period
          sin(2*pi*w/(365.25/2/7))*j_fac +
          cos(2*pi*w/(365.25/2/7))*j_fac +
          # adjustment for special weeks
          special_week*j_fac
      ),
      family = quasipoisson(link = 'log')
    ),
    'Serfling-Poisson GLM (single cycle)', 'serfl-1p', 'serfling', list(
      formula = formula(
        observed_deaths ~
          # log linear long term trend
          weeks_since_origin*j_fac +
          # seasonality
          # full year period
          sin(2*pi*w/(365.25/7))*j_fac +
          cos(2*pi*w/(365.25/7))*j_fac +
          # adjustment for special weeks
          special_week*j_fac +
          # exposures
          offset(log(exposure_pw_hmd))
      ),
      family = quasipoisson(link = 'log')
    ),
    'Poisson GAM', 'poisson-gam', 'gam', list(
      formula = formula(
        observed_deaths ~
          1 + sex + age_group +
          # log linear long term trend
          weeks_since_origin*j_fac +
          # penalized cyclic spline for seasonality
          s(w, bs = 'cp', k = 52, by = j_fac) +
          # adjustment for special weeks
          special_week*j_fac +
          # exposures
          offset(log(exposure_pw_hmd))
      ),
      family = quasipoisson(link = 'log'),
      method = 'REML'
    ),
    # 'Poisson GAM (no trend)', 'gam', list(
    #   formula = formula(
    #     observed_deaths ~
    #       1 + sex + age_group +
    #       # penalized cyclic spline for seasonality
    #       s(w, bs = 'cp', k = 52, by = j) +
    #       # adjustment for special weeks
    #       special_week*j +
    #       # exposures
    #       offset(log(exposure_pw_hmd))
    #   ),
    #   family = quasipoisson(link = 'log')
    # )
    'CODA LM', 'coda', 'coda', list(
    formula = formula(
      alr_share_on_annual_deaths ~
        # average alr proportion of deaths over the years
        # by epi-week, and stratified by age and sex
        w_fac*j_fac +
        # sex and age specific correlations between
        # centered and scaled "winter mortality" and
        # alr proportion of deaths by week
        mortality_winter_centered_yj:w_fac:j_fac
    ),
      reference_week = 0,
      winter_deaths_epi_weeks =
        IsoWeekToEpiWeek(1:9, w_start = 27),
      exogenous_epi_weeks =
        IsoWeekToEpiWeek(c(27:52, 1:9), w_start = 27)
    )
  ) %>%
  mutate(model_id = 1:n())

# Fit models and predict ------------------------------------------

# merge CV data with model specs
dat$fit_data <-
  dat$cv_data %>%
  nest(
    data = c(-country_code, -cv_id)
  ) %>%
  expand_grid(
    nest(
      dat$mod_spec,
      model_spec = c(-model_name, -model_name_short, -model_id, -model_class)
    )
  )

# This is a big loop which fits each model in the the model
# specifications to data for each country for each cross-validation
# series.
dat$fitted_models <-
  dat$fit_data %>%
  group_by(country_code, model_name, cv_id, model_name_short) %>%
  group_modify(.keep = FALSE, .f = ~{
    
    cat('Fit ', unlist(.y)[2],
        ' on CV set ', unlist(.y)[3],
        ' for ', unlist(.y)[1], '\n', sep = '')
    
    # prepare training and prediction data
    input_dat <- unnest(.x, data)
    dat_train <- filter(input_dat, sample == 'training')
    dat_pred <- input_dat
    
    # Average mortality models ----------------------------------------
    if (.x$model_class == 'avg') {
      
      # fit model and predict from it
      model_fit_and_predictions <-
        SimpleAveragesModel(
          df_training = dat_train, df_prediction = dat_pred,
          week_name = iso_week, death_name = observed_deaths,
          exposure_name = exposure_pw_hmd, year_name = year,
          n_years = .x$model_spec[[1]][[1]][[1]][['n_years']],
          stat = .x$model_spec[[1]][[1]][[1]][['stat']],
          sex, age_group
        )
      
    }
    
    # Serfling models -------------------------------------------------
    if (.x$model_class == 'serfling') {
      
      model_fit_and_predictions <-
        SerflingModel(
          df_training = dat_train, df_prediction = dat_pred,
          formula = .x$model_spec[[1]][[1]][[1]][['formula']],
          family = .x$model_spec[[1]][[1]][[1]][['family']]
        )
      
    }
    
    # GAM models ------------------------------------------------------
    if (.x$model_class == 'gam') {

      model_fit_and_predictions <-
        GAMModel(
          df_training = dat_train, df_prediction = dat_pred,
          formula = .x$model_spec[[1]][[1]][[1]][['formula']],
          family = .x$model_spec[[1]][[1]][[1]][['family']],
          method = .x$model_spec[[1]][[1]][[1]][['method']]
        )
            
    }

    # CODA models -----------------------------------------------------
    if (.x$model_class == 'coda') {
      
      model_fit_and_predictions <- CODAModel(
        df_input = input_dat,
        formula = .x$model_spec[[1]][[1]][[1]]$formula,
        winter_deaths_epi_weeks = .x$model_spec[[1]][[1]][[1]]$winter_deaths_epi_weeks,
        exogenous_epi_weeks = .x$model_spec[[1]][[1]][[1]]$exogenous_epi_weeks,
        reference_week = .x$model_spec[[1]][[1]][[1]]$reference_week,
        stratum_name = j, week_name = w, death_name = observed_deaths,
        exposure_name = exposure_pw_hmd, sample_name = sample, year_name = y
      ) 
        
    }
    
    # Assemble results ------------------------------------------------
    result <-
      tibble(
        fitted_model = list(model_fit_and_predictions$model),
        predictions = list(model_fit_and_predictions$predictions)
      )
    
    return(result)
    
  }) %>%
  ungroup()

# Plot observed vs. fitted ----------------------------------------

dat$fitted_models %>%
  select(-fitted_model) %>%
  group_by(country_code, model_name) %>%
  group_walk(~{
    
    fig_dat <-
      # predictions for single country and model
      .x %>%
      unnest(predictions) %>%
      group_by(cv_id, date, sample) %>%
      summarise(
        observed_deaths = sum(observed_deaths, na.rm = T),
        predicted_deaths = sum(predicted_deaths, na.rm = T)
      )
    
    scaler <- 2/max(fig_dat$observed_deaths)
    
    fig[[paste(.y[[1]], .y[[2]])]] <<-
      fig_dat %>%
      mutate(y_ob = cv_id + scaler*observed_deaths-1.5,
             y_pr = cv_id + scaler*predicted_deaths-1.5,
             y_to = scaler*observed_deaths-1.5) %>%
      ggplot(aes(x = date)) +
      geom_point(
        aes(
          group = cv_id,
          color = sample,
          y = y_ob
        ),
        size = 0.3
      ) +
      geom_point(aes(y = y_to), size = 0.3) +
      geom_line(aes(y = y_pr, group = cv_id, alpha = sample),
                color = 'red') +
      scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
      scale_y_continuous(breaks = 0:100, labels = c('Total', 1:100)) +
      scale_alpha_manual(values = c(training = 0.3, test = 1)) +
      guides(color = 'none', alpha = 'none') +
      scale_color_manual(values = glob$colors$sample) +
      glob$MyGGplotTheme() +
      labs(
        x = NULL, y = 'Cross Validation Series',
        title = paste(.y[[1]], .y[[2]])
      )
  })

# Exports ---------------------------------------------------------

fitted_models <- dat$fitted_models
save(
  fitted_models, file = paste0('out/', Sys.Date(), '-fitted_models.RData'),
  compress = 'xz'
)

ExportFiguresFromList(fig, path = 'out', add_date = TRUE)
