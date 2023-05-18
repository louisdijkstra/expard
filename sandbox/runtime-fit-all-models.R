library(expard)
library(tictoc)

tic()

cohort <- generate_cohort(
  n_patients = 1000,
  simulation_time = 100,
  n_drug_ADR_pairs = 1,
  risk_model = "risk_model_current_use()",
  min_chance_drug = 0.05,
  avg_duration = 5,
  prob_guaranteed_exposed = c(0),
  min_chance = 0.1,
  max_chance = .9,
  verbose = TRUE
)

pair <- cohort[[1]]

fit <- expard::fit_all_models(pair)
toc()