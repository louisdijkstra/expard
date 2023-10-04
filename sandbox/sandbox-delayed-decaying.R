library(expard)
library(dplyr)

logit <- function(p) { 
  log(p / (1 - p))  
}



# simulate some data -----------------------------------------------------------

cohort <- generate_cohort(
  n_patients = 1000,
  simulation_time = 100,
  n_drug_ADR_pairs = 1,
  risk_model = c("risk_model_delayed_decaying(2, 1, .5)"),
  min_chance_drug = .1,
  avg_duration = 5,
  max_chance_drug = NULL,
  prob_guaranteed_exposed = .1,
  min_chance = .1,
  max_chance = .5,
  verbose = TRUE
)

pair <- cohort[[1]]

determine_time_steps <- function(drug_history) { 
  simulation_time <- length(drug_history)
  sapply(1:simulation_time, function(t) { 
    # currently exposed or did not take the drug yet
    if (sum(drug_history[1:t]) == 0) { 
      return(0)  
    } else {
      time_steps_ago = t - min(which(drug_history[1:t] == 1))
      return(time_steps_ago) 
    }
  })
}


n_patients <- nrow(pair$drug_history)

time_steps <- do.call(rbind,
                      lapply(1:n_patients, function(i) {
                        determine_time_steps(pair$drug_history[i,])
                      }))

freq_table <- determine_frequency_unique_values(time_steps, pair$adr_history)


determine_loglikelihood_delayed_decaying <- function(param, 
                                             freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  mu_log <- param[3]
  sigma_log <- param[4]
  rate_log <- param[5]
  
  # transform to regular parameters
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  mu <- exp(mu_log)
  sigma <- exp(sigma_log)
  rate <- exp(rate_log)
  
  #combination <- delay(drug_history) + decay(drug_history)
  #combination / max(combination)
  
  # to make sure that the highest value is indeed 1
  normalizing_factor <- dnorm(mu, mu, sigma)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
    #risk_value = exp(-rate * (unique_value - 1))
    risk_value_delayed = dnorm(unique_value, mu, sigma) / normalizing_factor, 
    risk_value_decaying = exp(-rate * (unique_value - 1)), 
    risk_value = (risk_value_delayed + risk_value_decaying) / max(risk_value_delayed + risk_value_decaying)
  )
  
  # the risk when the time_since is 0 (never exposed or currently exposed)
  # is zero
  freq_table$risk_value[1] <- 0
  
  # determine the log-likelihood given the risks
  freq_table <- freq_table %>% rowwise() %>% 
    mutate(loglikelihood = -1*(n_adr*log((pi1 - pi0)*risk_value + pi0) + 
                                 n_no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))
  
  sum(freq_table$loglikelihood)
}

param <- c(-1, 0, log(1), log(1), log(1))



fit <- optim(param,
             determine_loglikelihood_delayed_decaying, 
             freq_table = freq_table,
             method = "Nelder-Mead",
             control = list(maxit = 1e8))

beta0 <- fit$par[1]
beta <- fit$par[2]
mu_log <- fit$par[3]
sigma_log <- fit$par[4]
rate_log <- fit$par[5]

exp(beta0) / (1 + exp(beta0))
exp(beta) / (1 + exp(beta)) 
exp(mu_log)
exp(sigma_log)
exp(rate_log)

f <- fit_all_models(pair, model = c("delayed", "decaying", "delayed+decaying"))

# 
# 
# 
# determine_frequency_unique_values <- function(mat, adr_history) {
#   unique_values <- sort(unique(as.vector(mat)))
#   
#   freq_table <-
#     lapply(sort(unique_values), function(unique_value) {
#       # where does this 'unique_value' occur in the data
#       indices <- which(mat == unique_value)
#       
#       # count how often an ADR occurred at those time points
#       n_adr <- sum(adr_history[indices])
#       
#       c(
#         unique_value = unique_value,
#         freq       = length(indices),
#         # how often does it occur
#         n_adr      = n_adr,
#         n_no_adr   = length(indices) - n_adr
#       )
#     })
#   
#   # create a tibble and return it
#   as_tibble(do.call(rbind, freq_table))
# }
# 
# 
# n_patients <- nrow(pair$drug_history)
# 
# 
# time_steps_delayed <- do.call(rbind,
#                       lapply(1:n_patients, function(i) {
#                         determine_time_steps_delayed(pair$drug_history[i,])
#                       }))
# 
# time_steps_decaying <- do.call(rbind,
#                               lapply(1:n_patients, function(i) {
#                                 determine_time_steps_decaying(pair$drug_history[i,])
#                               }))
# 
# r = sprintf("(%d,%d)",as.vector(time_steps_decaying), as.vector(time_steps_delayed))
# 
# unique(r)
# 
# View(r)
