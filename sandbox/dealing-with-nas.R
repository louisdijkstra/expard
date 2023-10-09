library(expard)

#' Apply Risk Model to Subset Drug History
#' 
#' Observations might not be present for certain patients, since they entered 
#' the data set later time = 1 or before the end. This function applies 
#' the risk model only to the observed time points. 
#' 
#' @param drug_history The drug history. Unobserved time points are denoted by 
#'                     \code{NA}
#' @param risk_model The risk model (Default \code{\link{risk_model_current_use}}). 
#' 
#' @return Vector of the same length as \code{drug_history}. None observed time points
#'         are denoted by \code{NA}
#' @export                 
apply_risk_model_only_to_observed_time <- function(drug_history, 
                                                   risk_model = expard::risk_model_current_use()) { 
  
  simulation_time <- length(drug_history)
  
  # determine the indices that are not NA
  indices_not_NA <- which(!is.na(drug_history))
  
  # get only the observed drug history given these indices
  observed_drug_history <- drug_history[indices_not_NA]

  # determine the risk for only this part of the drug history
  observed_risks <- risk_model(observed_drug_history)
  
  # initialize the risk vector with the same length as the original drug history
  risk <- rep(NA, simulation_time)
  
  # fill in the observed risks
  risk[indices_not_NA] <- observed_risks
  
  return(risk)
}

drug_history <- c(rep(NA, 5), rep(0,5), rep(1,5), rep(0,4), rep(NA,2))
drug_history <- c(rep(NA, 30), 0, rep(NA, 5))


apply_risk_model_only_to_observed_time(drug_history)

simulation_time <- length(drug_history)

indices_NA <- which(is.na(drug_history))
indices_not_NA <- setdiff(1:simulation_time, indices_NA)

observed_drug_history <- drug_history[indices_not_NA]

risk_model <- expard::risk_model_long_term(1,2)

observed_risks <- risk_model(observed_drug_history)
observed_risks

risk <- rep(NA, simulation_time)

risk[indices_not_NA] <- observed_risks


cohort <- generate_cohort(
  n_patients = 1000,
  simulation_time = 10,
  n_drug_ADR_pairs = 1,
  risk_model = c("risk_model_decaying(.9)"),
  min_chance_drug = .1,
  avg_duration = 5,
  max_chance_drug = NULL,
  prob_guaranteed_exposed = .1,
  min_chance = .1,
  max_chance = .8,
  verbose = TRUE
)

pair <- cohort[[1]]

pair$drug_history[10, 3:5] <- NA

pair$drug_history[12, 1:5] <- NA

pair$adr_history[10, 4:5] <- NA

pair$adr_history[12, 1:5] <- NA

expard::create2x2table(pair, method = "patient")

expard::create2x2table(pair, method = "time-point")

expard::create2x2table(pair, method = "drug-era")

pair

r <- fit_all_models(pair, models = c('no-association', 'current-use', 
                                     'past-use', 'withdrawal', 'delayed', 
                                     'decaying'))
