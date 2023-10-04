library(expard)
library(dplyr)

logit <- function(p) { 
  log(p / (1 - p))  
}



# simulate some data -----------------------------------------------------------

cohort <- generate_cohort(
  n_patients = 1000,
  simulation_time = 100,
  n_drug_ADR_pairs = 1,
  risk_model = c("risk_model_decaying_decaying(5, 2, .5)"),
  min_chance_drug = .1,
  avg_duration = 5,
  max_chance_drug = NULL,
  prob_guaranteed_exposed = .1,
  min_chance = .1,
  max_chance = .8,
  verbose = TRUE
)

pair <- cohort[[1]]




