#' Generate a Cohort 
#' 
#' Generates a cohort of patients, each with a 
#' drug prescription and ADR history. 
#' See for more information and details \code{\link{generate_patient}}.
#' 
#' @param n_patients Number of patients 
#' @inheritParams generate_patient
#' @param verbose Show progress bar (Default: \code{FALSE})
#' 
#' @return A \code{cohort} object; a list with 
#'       \item{\code{n_patients}}{The total number of patients}
#'       \item{\code{simulation_time}}{The total number of time steps}
#'       \item{\code{drug_history}}{A binary matrix with \code{n_patients}
#'             rows and \code{simulation_time} columns}
#'       \item{\code{adr_history}}{A binary matrix with \code{n_patients}
#'             rows and \code{simulation_time} columns}
#'       \item{\code{patient_profile}}{A list with the patient profiles for each patient}
#'
#' @seealso \code{\link{generate_patient}}   
#' @examples 
#' # choose a patient model (in this case based on the patients sex)
#' patient_model <- patient_model_sex(prob_male = 0.5)
#' (patient_profile <- patient_model())
#' 
#' # choose a risk model (effect of the drug is immediate)
#' risk_model <- risk_model_current_use()
#' 
#' # how the drug is prescribed over time
#' drug_model <- drug_model_markov_chain() 
#' 
#' # how the occurence of an ADR affect the prescription of the drug
#' adr_model = adr_model_no_effect()
#' 
#' # minimum and maximal probabilities for the drug being prescribed
#' min_chance_drug <- probability_model_sex(prob_male = .1, prob_female = .2)
#' max_chance_drug <- probability_model_sex(prob_male = .5, prob_female = .53)
#' 
#' # minimum and maximal probabilities for the ADR occurring 
#' min_chance_adr <- probability_model_sex(prob_male = .01, prob_female = .05)
#' max_chance_adr <- probability_model_sex(prob_male = .2, prob_female = .3)
#' 
#' # chance that the patient is exposed to the drug
#' prob_exposure <- probability_model_sex(prob_male = .9, prob_female = .7)
#' 
#' generate_cohort(n_patients = 10, 
#'                simulation_time = 20, 
#'                risk_model      = risk_model, 
#'                drug_model      = drug_model,
#'                adr_model       = adr_model, 
#'                min_chance_drug = min_chance_drug,
#'                max_chance_drug = max_chance_drug,
#'                min_chance_adr  = min_chance_adr, 
#'                max_chance_adr  = max_chance_adr, 
#'                patient_model   = patient_model,
#'                verbose = TRUE) 
#' @export
generate_cohort <- function(n_patients = 100, 
                            simulation_time = 100, 
                            n_drug_ADR_pairs = 50, 
                            risk_model = rep("risk_model_current_use()", n_drug_ADR_pairs), 
                            min_chance_drug = rep(.1, n_drug_ADR_pairs), 
                            avg_duration = rep(5, n_drug_ADR_pairs), 
                            max_chance_drug = rep(NULL, n_drug_ADR_pairs),
                            prob_guaranteed_exposed = rep(1, n_drug_ADR_pairs), 
                            min_chance  = rep(.1, n_drug_ADR_pairs), 
                            max_chance  = rep(.4, n_drug_ADR_pairs),
                            verbose = FALSE) { 
  
  
  if (verbose) { 
    cat("Generating the patients...\n")
    pb <- txtProgressBar(min = 0, max = n_patients, style = 3)
  }
  
  # generate patients 
  patients <- lapply(1:n_patients, function(i) { 
    
    guaranteed_exp <- rbinom(n_drug_ADR_pairs, 1, prob_guaranteed_exposed) 
    
    patient <- generate_patient(
      simulation_time,
      n_drug_ADR_pairs,
      risk_model,
      min_chance_drug,
      avg_duration,
      max_chance_drug,
      guaranteed_exposed = guaranteed_exp,
      min_chance,
      max_chance
    )
    
    if (verbose) { 
      setTxtProgressBar(pb, i)
    }
    return(patient)
  })
  
  if (verbose) { 
    close(pb) 
    cat("DONE generating the patients...\n\nOrganizing the data into matrices...\n")
    pb <- txtProgressBar(min = 0, max = n_drug_ADR_pairs, style = 3)
  }
  
  res <- lapply(1:n_drug_ADR_pairs, function(i) {
    drug_history <- matrix(rep(NA, simulation_time * n_patients), nrow = n_patients) 
    adr_history <- drug_history
    
    sapply(1:n_patients, function(p) {
      #print(patients[[p]][[i]])
      drug_history[p, ] <<- patients[[p]][[i]]$drug_history
      adr_history[p, ] <<- patients[[p]][[i]]$adr_history
    })

    if (verbose) { 
      setTxtProgressBar(pb, i)
    }
    
    list(drug_history = Matrix(drug_history, sparse = TRUE), 
         adr_history = Matrix(adr_history, sparse = TRUE))
  })

  if (verbose) { 
    close(pb) 
    cat("DONE organizing the data into matrices...\n")
  }
  
  res$n_patients <- n_patients
  res$n_drug_ADR_pairs <- n_drug_ADR_pairs
  res$simulation_time <- simulation_time
  
  class(res) <- "cohort"
  return(res)
}

#' Print function for \code{\link{generate_cohort}}
#' @export
print.cohort <- function(cohort) { 
  
  cat("Cohort\n\n")
  cat(sprintf("  No. patients:\t\t%d\n", cohort$n_patients))
  cat(sprintf("  No. drug-ADR-pairs:\t%d\n", cohort$n_drug_ADR_pairs))
  cat(sprintf("  No. time points:\t%d\n", cohort$simulation_time))
}