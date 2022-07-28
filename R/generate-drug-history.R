#' Generate Drug Prescription History
#' 
#' Returns whether or not the drug is prescribed 
#' at the next time point (\code{1}) or not (\code{0}). 
#' The probability depends on the drug prescriptions, 
#' the ADR history and the \code{drug_model}.
#' 
#' @param simulation_time Total number of time points
#' @param drug_model A drug model function, 
#'                           e.g., \code{\link{drug_model_markov_chain}}
#' @param min_chance The probability of the drug being prescribed
#'            when the drug history and the ADR history have no effect
#' @param max_chance The probability of the ADR when 
#'            the drug history and the ADR history have the highest 
#'            possible effect
#' 
#' @return \code{1} or \code{0}
#' @family Update functions
#' @examples 
#' # choose the drug prescription model
#' drug_model <- drug_model_markov_chain()
#' 
#' # choose probability models. Note that the probabilities change
#' # with sex
#' min_chance <- probability_model_sex(prob_male = .1, prob_female = .05)
#' max_chance <- probability_model_sex(prob_male = .7, prob_female = .8)
#' 
#' # create a patient profile: 
#' patient_model <- patient_model_sex(prob_male = 0.5) 
#' (patient_profile <- patient_model()) 
#' 
#' generate_drug_history(simulation_time = 10, 
#'                       drug_model, 
#'                       min_chance,
#'                       max_chance,
#'                       patient_profile)
#' @export
generate_drug_history <- function(simulation_time = 10, 
                                  drug_model,
                                  min_chance,
                                  max_chance, 
                                  patient_profile, ...) { 
  
  # determine first time the drug is prescribed 
  first <- determine_first_prescription(1, simulation_time, min_chance(patient_profile)) 
  
  if (first == simulation_time) { 
    drug_history <- c(rep(0, simulation_time - 1), 1) 
    return(drug_history)
  }
  
  # initialize the drug and ADR time processes 
  drug_history <- c(rep(0, first - 1), 1, rep(NA, simulation_time - first))
  
  
  # simulate the remaining time points
  for (t in (first+1):simulation_time) { 
    # get the probability of the next time point
    prob <- min_chance(patient_profile) + (max_chance(patient_profile) - min_chance(patient_profile)) * 
      drug_model(drug_history[1:(t-1)])
    
    drug_history[t] <- rbinom(1,1,prob)
  }
  
  return(drug_history)
}
