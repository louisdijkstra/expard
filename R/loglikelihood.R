#' Negative of the Log-Likelihood Function
#' 
#' \code{loglikelihood} returns the negative 
#' log-likelihood function, defined as 
#' \deqn{\pi = \frac{\exp(\beta_0 + \code{risks}\cdot \beta)}{1 + \exp(\beta_0 + \code{risks}\cdot \beta)}} 
#' where \eqn{\pi} is the probability. 
#' 
#' @param betas Vector with two entries: \code{beta0, beta}
#' @param risks Risks determined by the risk model. The risk 
#'              model was applied to the observed drug prescription 
#'              history
#' @export
loglikelihood <- function(betas,
                          risks,
                          adr_history) {

  # extract parameters
  beta0 <- betas[1]
  beta <- betas[2]

  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}


#' @export
loglikelihood_past <- function(param,
                               past,
                               drug_history,
                               adr_history) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  
  risk_model <- expard::risk_model_past(past)
  
  n_patients <- nrow(drug_history)
  simulation_time <- ncol(drug_history)
  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
  
  # go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(drug_history[i, ])
  }
  
  #risks <- Matrix(apply(drug_history, 1, function(d) risk_model(d)))
  
  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}


#' @export
loglikelihood_withdrawal <- function(param,
                            drug_history,
                            adr_history) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  phi <- param[3]
  
  risk_model <- expard::risk_model_withdrawal(exp(phi))
  
  
  n_patients <- nrow(drug_history)
  simulation_time <- ncol(drug_history)
  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
  
  # go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(drug_history[i, ])
  }
  
  #risks <- Matrix(apply(drug_history, 1, function(d) risk_model(d)))
  
  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}



#' @export
loglikelihood_delayed <- function(param,
                                  drug_history,
                                  adr_history) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  mu_log <- param[3]
  sigma_log <- param[4]
  
  risk_model <- expard::risk_model_delayed(exp(mu_log), exp(sigma_log))
  
  n_patients <- nrow(drug_history)
  simulation_time <- ncol(drug_history)
  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
  
  # go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(drug_history[i, ])
  }
  
  #risks <- Matrix(apply(drug_history, 1, function(d) risk_model(d)))
  
  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}





#' @export
loglikelihood_decaying <- function(param,
                                     drug_history,
                                     adr_history) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  phi <- param[3]
  
  risk_model <- expard::risk_model_decaying(exp(phi))
  
  
  n_patients <- nrow(drug_history)
  simulation_time <- ncol(drug_history)
  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
  
  # go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(drug_history[i, ])
  }
  
  #risks <- Matrix(apply(drug_history, 1, function(d) risk_model(d)))
  
  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}





#' @export
loglikelihood_delayed_decaying <- function(param,
                                  drug_history,
                                  adr_history) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  mu_log <- param[3]
  sigma_log <- param[4]
  rate_log <- param[5]
  
  risk_model <- expard::risk_model_delayed_decaying(exp(mu_log), exp(sigma_log), exp(rate_log))
  
  n_patients <- nrow(drug_history)
  simulation_time <- ncol(drug_history)
  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
  
  # go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(drug_history[i, ])
  }
  
  #risks <- Matrix(apply(drug_history, 1, function(d) risk_model(d)))
  
  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}


#risk_model_long_term <- function(rate, delay)

#' @export
loglikelihood_long_term <- function(param,
                                    drug_history,
                                    adr_history) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  rate_log <- param[3]
  delay_log <- param[4]
  
  risk_model <- expard::risk_model_long_term(exp(rate_log), exp(delay_log))
  
  n_patients <- nrow(drug_history)
  simulation_time <- ncol(drug_history)
  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
  
  # go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(drug_history[i, ])
  }
  
  #risks <- Matrix(apply(drug_history, 1, function(d) risk_model(d)))
  
  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}
