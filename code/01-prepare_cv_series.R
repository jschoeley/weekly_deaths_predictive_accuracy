# Prepare k-fold cross validation time-series splits of weekly
# deaths counts by age and sex
# Jonas Schöley

# Init ------------------------------------------------------------

source('code/00-init.R')

library(tidyverse)

dat <- list()
fig <- list()

# Constants -------------------------------------------------------

cnst <- list(
  # select population to work with
  country = 'DNK',
  # first week of test period
  iso_week_start_of_test = 10,
  # last week of test period
  iso_week_end_of_test = 33,
  # if mortality is elevated by a factor of <x> compared to prior week
  # or diminished by a factor of 1/<x> the week is marked as a "strange"
  # week
  strange_week_cutoff = 1.2
)

# Load data -------------------------------------------------------

# HMD short term mortality fluctuations data base
# https://www.mortality.org/Public/STMF/Outputs/stmf.csv
dat$raw <-
  readr::read_csv(
    "data/stmf.csv",
    col_types = "ciicddddddddddddlll",
    skip = 2, col_names = TRUE
  )

# Format data -----------------------------------------------------

dat$ts <-
  dat$raw %>%
  # select population of interest
  filter(
    CountryCode == cnst$country,
    # delete total sex category
    Sex != 'b'
  ) %>%
  # convert to long format over age
  select(
    year = Year, iso_week = Week,
    sex = Sex,
    D0_14:D85p, R0_14:R85p,
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
  rename(observed_deaths = D, observed_mortality = R) %>%
  # add variables relevant for analysis
  mutate(
    sex =
      factor(sex, names(glob$codebook$sex), glob$codebook$sex),
    age_group =
      factor(age_group, names(glob$codebook$age_group), glob$codebook$age_group),
    date =
      ISOWeekDate2Date(year, iso_week, 1),
    # date that corresponding epi-year started
    start_of_epi_year =
      ISOWeekDate2Date(
        ifelse(iso_week<=glob$week_epi_year_starts, year-1, year),
        glob$week_epi_year_starts, 1
      ),
    # current epi-year
    epi_year =
      paste0(lubridate::year(start_of_epi_year),'/',
             lubridate::year(start_of_epi_year)+1),
    # start of current epi-year
    epi_year_int =
      as.integer(substr(epi_year,1,4)),
    # weeks into epi-year (starting at 0)
    epi_week =
      difftime(date, start_of_epi_year, units = 'weeks') %>%
      floor() %>% as.integer(),
    # indicator variable for season
    season = case_when(
      iso_week %in% glob$seasons$northern$Winter ~ 'Winter',
      iso_week %in% glob$seasons$northern$Spring ~ 'Spring',
      iso_week %in% glob$seasons$northern$Summer ~ 'Summer',
      iso_week %in% glob$seasons$northern$Fall ~ 'Fall'
    ),
    # currently flu season?
    flu_season =
      iso_week %in% glob$flu_seasons$northern,
    # special weeks
    special_week =
      case_when(
        iso_week == 52 ~ 'Last week',
        iso_week == 1 ~ 'First week',
        TRUE ~ 'Regular week'
      ) %>%
      as_factor() %>%
      fct_relevel(
        c('Regular week', 'First week', 'Last week')
      )
  ) %>%
  # use same weekly exposure as used by HMD
  group_by(sex, age_group) %>%
  arrange(date) %>%
  mutate(
    exposure =
      observed_deaths / observed_mortality,
    # if mortality is 0 then exposures are NaN,
    # fill with previous value
    exposure =
      ifelse(observed_mortality == 0, NA, exposure)
  ) %>%
  fill(exposure, .direction = 'down') %>%
  ungroup()

# Check raw data --------------------------------------------------

# plot raw data
fig$raw_deaths <-
  dat$ts %>%
  ggplot() +
  geom_point(
    aes(x = date, y = observed_deaths),
    size = 0.2
  ) +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  scale_y_continuous(labels = scales::label_comma()) +
  facet_grid(age_group~sex, scales = 'free_y') +
  labs(y = 'Deaths per week', x = NULL) +
  glob$MyGGplotTheme(axis = 'xy')
fig$raw_exposures <-
  dat$ts %>%
  ggplot() +
  geom_point(
    aes(x = date, y = exposure),
    size = 0.2
  ) +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  scale_y_continuous(labels = scales::label_comma()) +
  facet_grid(age_group~sex, scales = 'free_y') +
  labs(y = 'Person-time exposure', x = NULL) +
  glob$MyGGplotTheme(axis = 'xy')

# Prepare cross-validation data sets ------------------------------

# cross-validation series
dat$cv <-
  # set up K-fold cross-validation
  #map2(2000+0:12, 2004+0:12, EpiYearSequence) # SWE
  #map2(2010+0:5, 2014+0:5, EpiYearSequence) # ENW
  map2(2006+0:9, 2010+0:9, EpiYearSequence) # DNK

# test-training data
dat$tt <-
  dat$cv %>%
  map(~filter(dat$ts, epi_year %in% .x)) %>%
  map(~mutate(.x, sample = ifelse(
    (year == max(epi_year_int)+1 & iso_week >= cnst$iso_week_start_of_test),
    'test', 'training')
  )) %>%
  bind_rows(.id = 'cv_id') %>%
  filter(
    sample == 'training' |
      (sample == 'test' & iso_week <= cnst$iso_week_end_of_test)
  ) %>%
  mutate(cv_id = as.integer(cv_id)) %>%
  arrange(cv_id, sex, age_group, date) %>%
  # add weeks since start of series
  group_by(cv_id) %>%
  mutate(
    weeks_since_origin =
      WeeksSinceOrigin(date, min(date))
  ) %>%
  ungroup() %>%
  select(
    # training test CV
    cv_id, sample,
    # strata
    sex, age_group,
    # time
    date, year, iso_week,
    epi_year, epi_year_int, start_of_epi_year, epi_week,
    # time flags
    season, flu_season, special_week,
    # population data
    observed_deaths, exposure, observed_mortality
  )

test_training_data <- dat$tt
save(test_training_data, file = 'out/test_training_data.RData')
write_csv(test_training_data, path = 'out/test_training_data.csv')

# Plot training-test split ----------------------------------------

# plot training-test split
fig$training_test_bars <-
  dat$tt %>%
  filter(age_group == '85+', sex == 'Male') %>%
  ggplot(aes(x = date, y = cv_id)) +
  geom_path(
    aes(color = sample, group = interaction(cv_id, sample)),
    size = 3
  ) +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  scale_y_continuous(breaks = 1:20) +
  labs(x = NULL, y = NULL) +
  guides(color = 'none') +
  scale_color_manual(values = glob$colors$sample) +
  glob$MyGGplotTheme()

# plot total weekly observed deaths
fig$weekly_observed_deaths <-
  dat$tt %>%
  group_by(date, cv_id) %>%
  summarise(
    observed_deaths = sum(observed_deaths)
  ) %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = observed_deaths)) +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  labs(x = NULL, y = 'Weekly observed deaths') +
  glob$MyGGplotTheme()

# plot K fold cross validation data series
fig$cross_validation <-
  dat$tt %>%
  group_by(date, cv_id, sample) %>%
  summarise(
    observed_deaths = sum(observed_deaths)
  ) %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = cv_id + 2/max(observed_deaths)*observed_deaths-1.5,
                 color = sample), size = 0.1) +
  geom_point(aes(y = 2/max(observed_deaths)*observed_deaths-1.5), size = 0.1) +
  glob$ggtheme +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  scale_y_continuous(breaks = 0:20, labels = c('Total', 1:20)) +
  labs(x = NULL, y = 'Cross Validation Series') +
  guides(color = 'none') +
  scale_color_manual(values = glob$colors$sample) +
  glob$MyGGplotTheme(grid = 'xy')

# Identify strange weeks ------------------------------------------

# identify strange weeks
fig$strange_weeks <-
  dat$ts %>%
  filter(year != 2020) %>%
  group_by(year, date, iso_week) %>%
  summarise(
    observed_deaths = sum(observed_deaths)
  ) %>%
  ungroup() %>%
  arrange(date) %>%
  mutate(
    diff_deaths = observed_deaths/lag(observed_deaths),
    strange =
      case_when(
        diff_deaths >= cnst$strange_week_cutoff ~ 'higher',
        diff_deaths <= 1/cnst$strange_week_cutoff ~ 'lower',
        TRUE ~ 'normal'
      )
  ) %>%
  filter(strange != 'normal') %>%
  ggplot(aes(x = iso_week, y = year)) +
  scale_y_continuous(breaks = 1900:2030) +
  scale_x_continuous(breaks = 1:52,
                     labels = function (x) ifelse(x%%5==0 | x == 1, x, '')) +
  geom_point(aes(fill = strange, color = strange), size = 5, shape = 21) +
  geom_point(aes(shape = strange, color = strange), size = 5,
             position = position_nudge(x = 0.02, y = 0.02)) +
  guides(color = 'none') +
  scale_shape_manual(values = c(higher = '+', lower = '-')) +
  scale_fill_manual(
    values = c(
      higher = glob$colors$discrete_light[2],
      lower = glob$colors$discrete_light[1]
    )
  ) +
  scale_color_manual(
    values = c(
      higher = glob$colors$discrete[2],
      lower = glob$colors$discrete[1]
    )
  ) +
  labs(x = 'Week into year', y = NULL) +
  glob$MyGGplotTheme(show_legend = FALSE, grid = 'xy')

# Export figures --------------------------------------------------

ExportFiguresFromList(fig, path = 'out', add_date = TRUE)
