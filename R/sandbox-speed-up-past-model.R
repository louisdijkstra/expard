library(expard)

library(microbenchmark)

n_patients = c(1000)
simulation_time = c(100)  
prob_exposed = c(0.5) 
avg_duration = c(5)
min_chance = c(1e-5) 
max_chance = c(.1)

min_chance_drug <- 1 - (1 - prob_exposed)^(1/simulation_time)

risk_model <- "current-use"

# first determine the risk model 
risk_model_fn <- switch(as.character(risk_model), 
                        "no-association" = "risk_model_no_association()", 
                        "current-use" = "risk_model_current_use()" , 
                        "past-use" = sprintf("risk_model_past(%d)", past) , 
                        "withdrawal" = sprintf("risk_model_withdrawal(%f)", rate) , 
                        "delayed" = sprintf("risk_model_delayed(%f, %f)", mu, sigma) , 
                        "decaying" = sprintf("risk_model_decaying(%f)", rate) , 
                        "delayed+decaying" = sprintf("risk_model_delayed_decaying(%f, %f, %f)", mu, sigma, rate) , 
                        "long-term" = sprintf("risk_model_long_term(%f, %f)", rate, delay))

# Generate a single drug-ADR pair
cohort <- generate_cohort(
  n_patients = n_patients,
  simulation_time = simulation_time,
  n_drug_ADR_pairs = 1,
  risk_model = risk_model_fn,
  min_chance_drug = min_chance_drug,
  avg_duration = avg_duration,
  prob_guaranteed_exposed = 0,
  min_chance = min_chance,
  max_chance = max_chance,
  verbose = FALSE
)

pair <- cohort[[1]]
sum(pair$drug_history)
sum(pair$adr_history)


past <- 5

drug_history <- pair$drug_history
adr_history <- pair$adr_history

risks <- drug_history




microbenchmark({
t(apply(drug_history, 1, function(drug_history_patient) {
  sapply(1:simulation_time, function(t) { 
    as.numeric( any(drug_history_patient[max(1,t-past):t] != 0))
  })
}))}, {# go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(drug_history[i, ])
  }}, times = 10)

risk_model <- risk_model_past(past)

risks <- matrix(0, nrow = n_patients, ncol = simulation_time)

# go over all patients 
for (i in 1:n_patients) { 
  # go over all timepoints 
  risks[i, ] <- risk_model(drug_history[i, ])
}

simulation_time <- length(drug_history) 

sapply(1:simulation_time, function(t) { 
  as.numeric( any(drug_history[max(1,t-past):t] != 0))
})


#risks <- Matrix(apply(drug_history, 1, function(d) risk_model(d)))

P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
-1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))


#expard::fit_all_models(pair)