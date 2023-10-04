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



determine_time_steps <- function(drug_history) { 
  simulation_time <- length(drug_history)
  
  # never took the drug
  if (sum(drug_history) == 0) { 
    return(rep(0, simulation_time))
  }
  
  sapply(1:simulation_time, function(t) { 
    # did not take the drug yet
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

exp(-rate * (0:10 - 1))

determine_loglikelihood_decaying <- function(param, 
                                            freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  rate_log <- param[3]
  
  # transform to regular parameters
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  rate <- exp(rate_log)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
    #risk_value = exp(-rate * (unique_value - 1))
    risk_value = exp(-rate * (unique_value - 1))
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

param <- c(-1, 0, log(1))



fit <- optim(param,
             determine_loglikelihood_decaying, 
             freq_table = freq_table,
             method = "Nelder-Mead",
             control = list(maxit = 1e8))

beta0 <- fit$par[1]
beta <- fit$par[2]
rate_log <- fit$par[3]

exp(beta0) / (1 + exp(beta0))
exp(beta) / (1 + exp(beta)) 
exp(rate_log)





fit <- optim(param,
             expard::loglikelihood_decaying, 
             freq_table = freq_table,
             method = "Nelder-Mead",
             control = list(maxit = 1e8))

beta0 <- fit$par[1]
beta <- fit$par[2]
rate_log <- fit$par[3]

exp(beta0) / (1 + exp(beta0))
exp(beta) / (1 + exp(beta)) 
exp(rate_log)

fit_model(pair, model = "decaying")





###########
# 
# determine_time_steps <- function(drug_history) { 
#   simulation_time <- length(drug_history)
#   
#   # never took the drug
#   if (sum(drug_history) == 0) { 
#     return(rep(0, simulation_time))
#   }
#   
#   sapply(1:simulation_time, function(t) { 
#     # did not take the drug yet
#     if (sum(drug_history[1:t]) == 0) { 
#       return(0)  
#     } else {
#       time_steps_ago = t - min(which(drug_history[1:t] == 1))
#       return(time_steps_ago) 
#     }
#   })
# }
# 
# n_patients <- nrow(pair$drug_history)
# 
# time_steps <- do.call(rbind,
#                       lapply(1:n_patients, function(i) {
#                         determine_time_steps(pair$drug_history[i,])
#                       }))
# 
# freq_table <- determine_frequency_unique_values(time_steps, pair$adr_history)
# 
# res <- optim(c(0,0,1),
#              expard::loglikelihood_decaying, 
#              freq_table = freq_table,
#              method = "Nelder-Mead",
#              control = list(maxit = maxiter))
# 
# beta0 <- res$par[1]
# beta <- res$par[2]
# rate_log <- res$par[3]
# 
# exp(beta0) / (1 + exp(beta0))
# exp(beta) / (1 + exp(beta)) 
# exp(rate_log)
# 
# 
fit_model(pair, model = "decaying")

f = fit_all_models(pair, models = c("no-association", "current-use", "withdrawal", "delayed", "decaying"))
