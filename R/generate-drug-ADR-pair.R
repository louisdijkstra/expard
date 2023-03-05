#' Generate a Drug-ADR pair
#' 
#' Generates the drug- and ADR history
#' 
#' @param simulation_time The total number of time steps
#' @param risk_model One of the risk models
#' @param min_chance_drug The probability of the drug being prescribed
#'            when the drug history and the ADR history have no effect
#' @param avg_duration Average number of time points a patient is exposed
#'                     once exposed to the drug (Determines \code{max_chance})
#' @param max_chance_drug The probability of the ADR when 
#'            the drug history and the ADR history have the highest 
#'            possible effect (Default: \code{NULL})
#' @param guaranteed_exposed If \code{TRUE}, the patient is exposed to the drug
#'            at least once
#' @param min_chance The probability of the ADR when 
#'            the drug history has no effect 
#' @param max_chance The probability of the ADR when 
#'            the drug history has the highest possible effect
#' 
#' @return A list with
#'       \item{\code{drug_history}}{binary vector of the 
#'                drug history}
#'       \item{\code{adr_history}}{binary vector of the 
#'                ADR history}
#'
#' @seealso \code{\link{generate_cohort}}  
#' @examples 
#' @export
generate_drug_ADR_pair <- function(simulation_time = 100, 
                             risk_model = risk_model_current_use(), 
                             min_chance_drug = .1, 
                             avg_duration = 5, 
                             max_chance_drug = NULL,
                             guaranteed_exposed = TRUE, 
                             min_chance  = .1, 
                             max_chance  = .4) { 
  
  
  drug_history <- generate_drug_history(simulation_time, min_chance_drug, 
                                        avg_duration, max_chance_drug, 
                                        guaranteed_exposed)
  
  
  adr_history <- generate_adr_history(drug_history, risk_model, min_chance, max_chance) 
  
  res <- list(
    drug_history = drug_history,
    adr_history = adr_history
  )
  return(res)
}

#' Function for printing a patient 
#' @export
print.patient <- function(patient) { 
  cat(sprintf("Patient\n")) 
  if (length(patient$patient_profile) != 0) {
    for (i in 1:length(patient$patient_profile)) { 
      cat(sprintf("\t-- %s: ", names(patient$patient_profile)[i]))
      cat(patient$patient_profile[[i]])
      cat(sprintf("\n"))
    }
  }
  cat(sprintf("\nNo. of time points: %d\n\n", patient$simulation_time))
  
  cat(sprintf("drug: "))
  for(i in 1:patient$simulation_time) { 
    if (patient$drug_history[i] == 1) { 
      cat(green(1))
    } else { 
      cat(blue("."))
    }
  }
  cat(sprintf("\nADR:  "))
  for(i in 1:patient$simulation_time) { 
    if (patient$adr_history[i] == 1) { 
      cat(red(1))
    } else { 
      cat(blue("."))
    }
  }
  cat(sprintf("\n"))
}
