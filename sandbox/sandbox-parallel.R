library(parallel)









library(expard)
library(dplyr)


# simulate some data -----------------------------------------------------------



cohort <- generate_cohort(
  n_patients = 5000,
  simulation_time = 56,
  n_drug_ADR_pairs = 1,
  risk_model = c("risk_model_withdrawal(.5)"),
  min_chance_drug = .1,
  avg_duration = 5,
  max_chance_drug = NULL,
  prob_guaranteed_exposed = .1,
  min_chance = .01,
  max_chance = .6,
  verbose = TRUE
)

pair <- cohort[[1]]


library(tictoc) 

tic()
r_par <- expard::fit_all_models(pair, models = models, mc.cores = 15)
toc()

tic()
r_seq <- expard::fit_all_models(pair, models = models, mc.cores = 1)
toc()

models = c(
  'no-association',
  'current-use',
  'past-use',
  'withdrawal',
  'delayed',
  'decaying',
  'delayed+decaying',
  'long-term'
)

expard::fit_model(pair, model = "past-use(20)")
toc()
