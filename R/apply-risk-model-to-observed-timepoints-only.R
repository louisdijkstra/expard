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