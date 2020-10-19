# Short term mortality forecast model evaluation
Jonas Schöley

We evaluate the ability of various models to predict weekly death counts into the near future. The basis for the evaluation are 5 overlapping time series of 7 years training data, followed by less than a year of test data. Data for multiple countries are available.

- specify models to test in `code/00-models.R`
- fit the models via `code/01-fit_models.R`
- generate a report by running `code/predictive_performance_report.Rmd`

Cross-validation data is found in `data/2020-10-05-xstmf_cv.RData` and derived from the HMD's short term mortality fluctuations data base at https://www.mortality.org/Public/STMF/Outputs/stmf.csv.