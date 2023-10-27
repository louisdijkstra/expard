library(expard)
library(dplyr)

logit <- function(p) { 
  log(p / (1 - p))  
}


# simulate some data -----------------------------------------------------------

rate <- .5

cohort <- generate_cohort(
  n_patients = 1000,
  simulation_time = 100,
  n_drug_ADR_pairs = 1,
  risk_model = c(sprintf("risk_model_withdrawal(%g)", rate)),
  min_chance_drug = .1,
  avg_duration = 5,
  max_chance_drug = NULL,
  prob_guaranteed_exposed = .1,
  min_chance = .1,
  max_chance = .8,
  verbose = TRUE
)

pair <- cohort[[1]]

# process the data -------------------------------------------------------------
# determine times since last exposure for a drug history of a single patient
determine_time_steps_ago <- function(drug_history) {
  sapply(1:length(drug_history), function(t) { 
    # currently exposed or did not take the drug yet
    if (drug_history[t] == 1 || sum(drug_history[1:t]) == 0) { 
      return(0)  
    } else {
      time_steps_ago = t - max(which(drug_history[1:t] == 1))
      return(time_steps_ago) 
    }
  })
}

#apply_function_to_observed_timepoints(pair$drug_history[10,], 
#                                      determine_time_steps_ago)

n_patients <- nrow(pair$drug_history)

time_steps_ago <- do.call(rbind, 
                          lapply(1:n_patients, function(i) { 
                            apply_function_to_observed_timepoints(pair$drug_history[i, ], 
                                                                  determine_time_steps_ago)
                            #determine_time_steps_ago(pair$drug_history[i, ])
                          }))

freq_table <- expard::determine_frequency_unique_values(time_steps_ago, pair$adr_history)
# determine times since last exposure for a drug history of a single patient
determine_time_steps_ago <- function(drug_history) {
  sapply(1:length(drug_history), function(t) { 
    # currently exposed or did not take the drug yet
    if (drug_history[t] == 1 || sum(drug_history[1:t]) == 0) { 
      return(0)  
    } else {
      time_steps_ago = t - max(which(drug_history[1:t] == 1))
      return(time_steps_ago) 
    }
  })
}

#apply_function_to_observed_timepoints(pair$drug_history[10,], 
#                                      determine_time_steps_ago)

n_patients <- nrow(pair$drug_history)

time_steps_ago <- do.call(rbind, 
                          lapply(1:n_patients, function(i) { 
                            apply_function_to_observed_timepoints(pair$drug_history[i, ], 
                                                                  determine_time_steps_ago)
                            #determine_time_steps_ago(pair$drug_history[i, ])
                          }))

freq_table <- expard::determine_frequency_unique_values(time_steps_ago, pair$adr_history)



freq_table <- freq_table %>% mutate(
  risk_value = exp(-rate * (unique_value - 1))
)

freq_table$risk_value[1] <- 0

freq_table$r <- freq_table$risk_value / sum(freq_table$risk_value)

fit_model(pair, model = "withdrawal")

sum((f$risk_value) * f$n_adr / f$s + (1 - f$risk_value) * f$n_no_adr / f$s)



sum(freq_table$n_adr / (freq_table$n_adr + freq_table$n_no_adr))

f <- freq_table
f$s <- f$n_adr + f$n_no_adr

sum(f$n_adr / (f$s) * f$risk_value + (1 - f$n_no_adr / (f$s)))

sum(freq_table$n_no_adr / (freq_table$n_adr + freq_table$n_no_adr) * freq_table$risk_value)


freq_table$n_adr / (freq_table$n_adr + freq_table$n_no_adr) * freq_table$r


