library(expard)


model <- c("risk_model_past(5)")
maxiter = 100
SIMULATION_TIME <- 100
N_PATIENTS <- 100

cohort <- generate_cohort(
  n_patients = N_PATIENTS,
  simulation_time = SIMULATION_TIME,
  n_drug_ADR_pairs = 1,
  risk_model = model,
  min_chance_drug = .01,
  avg_duration = 5,
  max_chance_drug = NULL,
  prob_guaranteed_exposed = .1,
  min_chance = .1,
  max_chance = .4,
  verbose = TRUE
)
pair <- cohort[[1]]

# initialize the fit --------------------------------

# past <- 1:(fit$simulation_time - 1)
past <- c(1,3,5,7)

refit_using_old_version <- TRUE

if (refit_using_old_version) { 

fit <- data.frame(
  expand.grid(
    n_patients = nrow(pair$drug_history),
    simulation_time = ncol(pair$drug_history),
    model = model[1],
    n_param = 3,
    loglikelihood = NA,
    converged = NA, 
    past = past
  )
)

cat(sprintf("Parameter 'past': \n"))
pb <- txtProgressBar(min = 1, max = max(past), style = 3)  

estimates <- lapply(past, function(d) { 
  res <- optim(c(0,0,-1),
               loglikelihood_past, 
               past = d, 
               drug_history = pair$drug_history,
               adr_history = pair$adr_history,
               method = "Nelder-Mead",
               control = list(maxit = maxiter))
  setTxtProgressBar(pb, d)
  return(res)
})

close(pb)

fit$loglikelihood <- sapply(estimates, function(est) est$value)
fit$p0 <- sapply(estimates, function(est) {
  beta0 <- est$par[1]
  exp(beta0) / (1 + exp(beta0))
})
fit$p1 <- sapply(estimates, function(est) {
  beta <- est$par[2]
  exp(est$par[1] + beta) / (1 + exp(est$par[1] + beta))
})
fit$converged <- sapply(estimates, function(est) est$convergence == 0)

fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
fit$bestBIC <- min(fit$BIC)
}

## New approach 

n_patients <- nrow(pair$drug_history)
simulation_time <- ncol(pair$drug_history)

#index_patients_not_exposed <- which(rowSums(pair$drug_history) == 0)
#n0. <- simulation_time * length(index_patients_not_exposed)
#index_patients_exposed <- setdiff(1:n_patients, index_patients_not_exposed) 

#n1. <- simulation_time * length(index_patients_exposed)
#n.1 <- sum(pair$adr_history)
#n_ADR_not_exposed_patients <- sum(pair$adr_history[index_patients_not_exposed, ])

#n01 <- n_ADR_not_exposed_patients


past <- 1





########################################################

# select the risk model 
risk_model <- expard::risk_model_past(past)

# 'convert' the drug prescriptions. They reflect which period is considered
# to have an increased risk
risks <- matrix(0, nrow = n_patients, ncol = simulation_time)

# go over all patients 
for (i in 1:n_patients) { 
  # go over all timepoints 
  risks[i, ] <- risk_model(pair$drug_history[i, ])
}

#' given the risk, determine the 2x2 table: 
#' 
#'                  ADR
#'         |     1       0      |        
#'       ------------------------------------
#' risk  1 |    n11     n10     | n1.
#'       0 |    n01     n00     | n0. 
#'       ------------------------------------
#'         |    n.1     n.0     | simulation_time*n_patients

# determine the marginals
n1. <- sum(risks)
n0. <- simulation_time*n_patients - n1. 
n.1 <- sum(pair$adr_history)
n.0 <- simulation_time*n_patients - n.1 

# determine the entries 
n11 <- sum(risks & pair$adr_history)
n01 <- n.1 - n11
n10 <- n1. - n11
n00 <- n0. - n01

# determine the estimates
pi1 <- n11 / n1. 
pi0 <- n01 / n0. 

fit$n_param <- 3
fit$loglikelihood <- -1*n11*log(pi1) - n10*log(1 - pi1) - n01*log(pi0) - n00*log(1 - pi0)
fit$converged <- TRUE

fit$p1 <- pi1
fit$p0 <- pi0

fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
fit$bestBIC <- min(fit$BIC)

return(fit)










past <- 1:(simulation_time - 1)

fit <- data.frame(
  expand.grid(
    n_patients = nrow(pair$drug_history),
    simulation_time = ncol(pair$drug_history),
    model = model[1],
    n_param = 3,
    loglikelihood = NA,
    converged = NA, 
    past = past
  )
)

cat(sprintf("Parameter 'past': \n"))
pb <- txtProgressBar(min = 1, max = max(past), style = 3)  

estimates <- lapply(past, function(d) { 
  # select the risk model 
  risk_model <- expard::risk_model_past(d)
  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
  
  # go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(pair$drug_history[i, ])
  }
  
  #' given the risk, determine the 2x2 table: 
  #' 
  #'                  ADR
  #'         |     1       0      |        
  #'       ------------------------------------
  #' risk  1 |    n11     n10     | n1.
  #'       0 |    n01     n00     | n0. 
  #'       ------------------------------------
  #'         |    n.1     n.0     | simulation_time*n_patients
  
  # determine the marginals
  n1. <- sum(risks)
  n0. <- simulation_time*n_patients - n1. 
  n.1 <- sum(pair$adr_history)
  n.0 <- simulation_time*n_patients - n.1 
  
  # determine the entries 
  n11 <- sum(risks & pair$adr_history)
  n01 <- n.1 - n11
  n10 <- n1. - n11
  n00 <- n0. - n01
  
  setTxtProgressBar(pb, d)
  return(list(p1 = n11 / n1., 
              p0 = n01 / n0., 
              value = -1*n11*log(pi1) - n10*log(1 - pi1) - n01*log(pi0) - n00*log(1 - pi0), 
              converged = TRUE))
})

close(pb)

fit$loglikelihood <- sapply(estimates, function(est) est$value)
fit$p0 <- sapply(estimates, function(est) est$p0 )
fit$p1 <- sapply(estimates, function(est) est$p1 )
fit$converged <- sapply(estimates, function(est) est$convergence == 0)

fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
fit$bestBIC <- min(fit$BIC)



# 
# 
# #n.0 <- simulation_time * n_patients - n.1 
# #n0. <- simulation_time * n_patients - n1. 
# n10 <- n1. - n11
# n00 <- n.0 - n10
# 
# pi1 <- n11 / n1. 
# pi0 <- n01 / n0. 
# cat(sprintf("pi0: %g \t pi1: %g\n", pi0, pi1))
# 
# 
# #logl <- length(index_patients_not_exposed) * SIMULATION_TIME * log(1 - )
