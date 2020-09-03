# Model weekly death counts by sex and age
# Jonas Schöley

# Init ------------------------------------------------------------

source('code/00-init.R')

library(tidyverse)
library(mgcv)

dat <- list()
fig <- list()

# Load data -------------------------------------------------------

load('out/test_training_data.RData')

# Specifications of models to test --------------------------------

# just a big list specifying all the models to be tested

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
          sin(2*pi*epi_week/(365.25/7))*sex*age_group +
          cos(2*pi*epi_week/(365.25/7))*sex*age_group +
          # half year period
          sin(2*pi*epi_week/(365.25/2/7))*sex*age_group +
          cos(2*pi*epi_week/(365.25/2/7))*sex*age_group +
          # adjustment for special weeks
          special_week*sex*age_group +
          # exposures
          offset(log(exposure))
      ),
      family = quasipoisson(link = 'log')
    ),
    'Serfling-Poisson GLM\n(only full year cycle)', 'glm', list(
      formula = formula(
        observed_deaths ~
          # log linear long term trend
          weeks_since_origin*sex*age_group +
          # seasonality
          # full year period
          sin(2*pi*epi_week/(365.25/7))*sex*age_group +
          cos(2*pi*epi_week/(365.25/7))*sex*age_group +
          # adjustment for special weeks
          special_week*sex*age_group +
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
          s(epi_week, bs = 'cp', k = 52, by = sex_age_interaction) +
          # adjustment for special weeks
          special_week*sex*age_group +
          # exposures
          offset(log(exposure))
      ),
      family = quasipoisson(link = 'log')
    ),
    'Poisson GAM (no trend)', 'gam', list(
      formula = formula(
        observed_deaths ~
          1 + sex + age_group +
          # penalized cyclic spline for seasonality
          s(epi_week, bs = 'cp', k = 52, by = sex_age_interaction) +
          # adjustment for special weeks
          special_week*sex*age_group +
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
          s(epi_week, bs = 'cp', k = 52, by = sex_age_interaction) +
          # adjustment for special weeks
          special_week*sex*age_group +
          # exposures
          offset(log(exposure))
      ),
      family = nb(link = 'log')
    ),
    'Neg-Bin. GAM-RE', 'gam', list(
      formula = formula(
        observed_deaths ~
          1 + sex + age_group +
          # log linear long term trend
          weeks_since_origin*sex_age_interaction +
          # penalized cyclic spline for seasonality
          s(epi_week, bs = 'cp', k = 52, by = sex_age_interaction) +
          # adjustment for special weeks
          special_week*sex*age_group +
          # flu season year random effect adjustment
          s(flu_season_year_sex_age_interaction, bs = 're') +
          # exposures
          offset(log(exposure))
      ),
      family = nb(link = 'log')
    )
  ) %>%
  mutate(model_id = 1:n())

# Fit models and predict ------------------------------------------

dat$tt <-
  test_training_data %>%
  nest(
    data = c(-cv_id)
  ) %>%
  expand_grid(
    nest(
      dat$mod_spec,
      model_spec = c(-model_name, -model_id, -model_class)
    )
  )

dat$mod <-
  dat$tt %>%
  group_by(model_name, cv_id) %>%
  group_modify(~{
    
    cat('Fit ', unlist(.y)[1], ' on CV set ', unlist(.y)[2], '\n', sep = '')
    
    ### PREPARE INPUT DATA
    
    # season-years in training data
    # I use this later to check if prediction
    # outside the range of epi_years in training
    # data will be required
    epi_years_in_training <-
      .x$data[[1]] %>%
      filter(sample == 'training') %>%
      pull(epi_year) %>% unique()
    
    # prepare for model
    model_data <-
      .x$data[[1]] %>%
      mutate(
        sex_age_interaction =
          interaction(sex, age_group),
        include_re =
          ifelse(epi_year %in% epi_years_in_training, 1, 0),
        flu_season_year_sex_age_interaction =
          interaction(flu_season, epi_year, sex, age_group)
      )
    
    training_data <-
      model_data %>%
      filter(sample == 'training')
    
    ### TRAIN MODELS
    
    if (.x$model_class == 'avgmx') {
      
      # calculate average weekly mortality rate over training data
      model <-
        AverageMortalityModel(
          training_data,
          iso_week,
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
          newdata$epi_year_sex_age_interaction <- NULL
          exclude <- 's(epi_year_sex_age_interaction)'
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
            predicted_mortality * exposure
        )
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

mod <- dat$mod
save(mod, file = './out/model_fits.RData', compress = 'xz')

# Plot observed vs. fitted ----------------------------------------

dat$mod_spec %>%
  group_by(model_name) %>%
  group_walk(~{
    
    fig_dat <-
      dat$mod %>%
      filter(model_name == .y[[1]]) %>%
      unnest(predictions) %>%
      group_by(model_name, cv_id, date, sample) %>%
      summarise(
        observed_deaths = sum(observed_deaths, na.rm = T),
        predicted_deaths = sum(predicted_deaths, na.rm = T)
      )
    
    scaler <- 2/max(fig_dat$observed_deaths)
    
    fig[[.y[[1]]]] <<-
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
      facet_wrap(~model_name, nrow = 2) +
      scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
      scale_y_continuous(breaks = 0:100, labels = c('Total', 1:100)) +
      scale_alpha_manual(values = c(training = 0.3, test = 1)) +
      labs(x = NULL, y = 'Cross Validation Series') +
      guides(color = 'none', alpha = 'none') +
      scale_color_manual(values = glob$colors$sample) +
      glob$MyGGplotTheme()
  })

ExportFiguresFromList(fig, path = 'out', add_date = TRUE)
