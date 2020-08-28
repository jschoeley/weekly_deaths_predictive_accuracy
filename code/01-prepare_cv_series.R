# Prepare k-fold cross validation time-series splits of weekly
# deaths counts by age and sex
# Jonas Schöley

# Init ------------------------------------------------------------

source('code/00-init.R')

library(tidyverse)

dat <- list()
fig <- list()

# Load data -------------------------------------------------------

# HMD short term mortality fluctuations data base
# https://www.mortality.org/Public/STMF/Outputs/stmf.csv
dat$raw <-
  readr::read_csv(
    "data/stmf-2.csv",
    col_types = "ciicddddddddddddlll",
    skip = 2, col_names = TRUE
  )

# Format data -----------------------------------------------------

dat$ts <-
  dat$raw %>%
  filter(
    CountryCode == glob$country,
    # exclude total sex category
    Sex != 'b'
  ) %>%
  select(
    year = Year, iso_week = Week, sex = Sex, D0_14:D85p, R0_14:R85p,
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
  mutate(
    sex =
      factor(sex, names(glob$codebook$sex), glob$codebook$sex),
    age_group =
      factor(age_group, names(glob$codebook$age_group), glob$codebook$age_group),
    # use same weekly exposure as used by HMD
    exposure =
      observed_deaths / observed_mortality,
    # assuming that HMD uses iso-weeks
    date =
      ISOWeekDate2Date(year, iso_week, 1),
    start_of_epi_year =
      ISOWeekDate2Date(
        ifelse(iso_week<=glob$week_epi_year_starts, year-1, year),
        glob$week_epi_year_starts, 1
      ),
    epi_year =
      paste0(lubridate::year(start_of_epi_year),'/',
             lubridate::year(start_of_epi_year)+1),
    epi_year_int =
      as.integer(substr(epi_year,1,4)),
    epi_week =
      difftime(date, start_of_epi_year, units = 'weeks') %>%
      floor() %>% as.integer(),
    season = case_when(
      iso_week %in% glob$seasons$northern$Winter ~ 'Winter',
      iso_week %in% glob$seasons$northern$Spring ~ 'Spring',
      iso_week %in% glob$seasons$northern$Summer ~ 'Summer',
      iso_week %in% glob$seasons$northern$Fall ~ 'Fall'
    ),
    flu_season =
      iso_week %in% glob$flu_seasons$northern,
    last_iso_week =
      ifelse(iso_week == 52, 1, 0),
    first_iso_week =
      ifelse(iso_week == 1, 1, 0)
  )

# Check raw data --------------------------------------------------

# plot raw data
fig$raw_data <-
  dat$ts %>%
  ggplot() +
  geom_point(
    aes(x = date, y = observed_deaths/exposure,
        color = age_group),
    size = 0.2
  ) +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  scale_y_log10() +
  facet_wrap(~sex, nrow = 2) +
  glob$ggtheme

# Prepare cross-validation data sets ------------------------------

# cross-validation series
dat$cv <-
  # set up K-fold cross-validation
  #map2(2000+0:12, 2004+0:12, EpiYearSequence) # SWE
  map2(2010+0:5, 2014+0:5, EpiYearSequence) # ENW
  #map2(2006+0:9, 2010+0:9, EpiYearSequence) # DNK

# test-training data
dat$tt <-
  dat$cv %>%
  map(~filter(dat$ts, epi_year %in% .x)) %>%
  map(~mutate(
    .x, training = ifelse(
        (year == max(epi_year_int)+1 & iso_week >= 10),
      'test', 'training')
  )) %>%
  bind_rows(.id = 'cv_id') %>%
  mutate(cv_id = as.integer(cv_id))

test_training_data <- dat$tt
save(test_training_data, file = 'out/test_training_data.RData')

# Plot training-test split ----------------------------------------

# plot training-test split
fig$training_test_bars <-
  dat$tt %>%
  filter(age_group == '85+', sex == 'Male') %>%
  ggplot(aes(x = date, y = cv_id)) +
  geom_path(
    aes(color = training, group = interaction(cv_id, training)),
    size = 3
  ) +
  glob$ggtheme +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  scale_y_continuous(breaks = 1:20) +
  labs(x = NULL, y = NULL) +
  guides(color = 'none') +
  scale_color_manual(values = glob$colors$col_training)

# plot total weekly observed deaths
fig$weekly_observed_deaths <-
  dat$tt %>%
  group_by(date, cv_id) %>%
  summarise(
    observed_deaths = sum(observed_deaths)
  ) %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = observed_deaths)) +
  glob$ggtheme +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  labs(x = NULL, y = 'Weekly observed deaths')

# plot K fold cross validation data series
fig$cross_validation <-
  dat$tt %>%
  group_by(date, cv_id, training) %>%
  summarise(
    observed_deaths = sum(observed_deaths)
  ) %>%
  ggplot(aes(x = date)) +
  geom_point(aes(y = cv_id + 0.0001*observed_deaths-1,
                 color = training), size = 0.1) +
  geom_point(aes(y = 0.0001*observed_deaths-1), size = 0.1) +
  glob$ggtheme +
  scale_x_date(date_breaks = '1 year', date_labels = '%Y') +
  scale_y_continuous(breaks = 0:20, labels = c('Total', 1:20)) +
  labs(x = NULL, y = 'Cross Validation Series') +
  guides(color = 'none') +
  scale_color_manual(values = glob$colors$col_training) +
  theme(
    panel.grid.major.y = element_line(color = 'grey80', linetype = 1, size = 0.2),
    panel.grid.major.x = element_line(color = 'grey80', linetype = 1, size = 0.2)
  )

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
      ifelse(diff_deaths >= 1.25 | diff_deaths <= 0.8, TRUE, FALSE)
  ) %>%
  ggplot(aes(x = iso_week, y = diff_deaths)) +
  geom_vline(xintercept = 1:52, size = 0.2, color = 'grey80') +
  geom_hline(yintercept = 1, color = 'grey80') +
  geom_point(aes(color = strange)) +
  guides(color = 'none') +
  scale_y_log10(
    breaks = c(0.8, 1, 1.25), labels = c('100:125', '0', '125:100')
  ) +
  facet_wrap(~year) +
  glob$ggtheme

# Export figures --------------------------------------------------

ExportFiguresFromList(fig, path = 'out', add_date = TRUE)
