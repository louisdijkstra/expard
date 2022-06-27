#' Check Cohort Object
#' 
#' \code{check_cohort} checks whether the argument for a cohort 
#' is valid and turns the object into a \code{cohort} object from 
#' the \code{hccdsim} package, which allows us to use 
#' the output function. In case the input is invalid, 
#' the code halts with an appropriate message. 
#' The minimum requirement is that the \code{cohort}
#' is a list with the items: \code{drug_history}
#' and \code{adr_history}. Both should be 
#' \emph{binary} matrix with \code{n_patients} rows and 
#' \code{simulation_time} columns.  
#' 
#' @param cohort A cohort 
#' 
#' @return A \code{cohort} object, see \code{\link{generate_cohort}}  
#'      
#' @seealso \code{\link{generate_cohort}}
#' @export
check_cohort <- function(cohort) { 
  
  if (!("drug_history" %in% names(cohort)) ) { 
    stop("cohort should contain drug_history")
  }
  
  if (!("adr_history" %in% names(cohort)) ) { 
    stop("cohort should contain adr_history")
  }
  
  # check whether drug_history and adr_history
  # are matrices
  dim_drug_history <- dim(cohort$drug_history)
  if (is.null(dim_drug_history)) { 
    stop("drug_history should be a matrix") 
  }
  
  dim_adr_history <- dim(cohort$adr_history) 
  if (is.null(dim_adr_history)) { 
    stop("adr_history should be a matrix") 
  }
  
  # check whether drug_history and adr_history
  # have the same dimensions
  if (dim_drug_history[1] != dim_adr_history[1]) { 
    stop(sprintf("number of patients differ between drug_history 
       and adr_history, (%d and %d, respectively)", 
                 dim_drug_history[1], 
                 dim_adr_history[1])) 
  }
  
  if (dim_drug_history[2] != dim_adr_history[2]) { 
    stop(sprintf("number of time points differ between drug_history 
       and adr_history, (%d and %d, respectively)", 
                 dim_drug_history[2], 
                 dim_adr_history[2])) 
  }
  
  # turn the cohort object into a cohort object 
  cohort$n_patients <- dim_drug_history[1]
  cohort$simulation_time <- dim_drug_history[2]
  
  class(cohort) <- "cohort"
  return(cohort)
}