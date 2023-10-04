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
  risk_model = c("risk_model_delayed(5, 2)"),
  min_chance_drug = .1,
  avg_duration = 5,
  max_chance_drug = NULL,
  prob_guaranteed_exposed = .1,
  min_chance = .01,
  max_chance = .6,
  verbose = TRUE
)

pair <- cohort[[1]]



# process the data -------------------------------------------------------------

determine_time_steps_since_start <- function(drug_history) { 
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


time_steps_since_start <- do.call(rbind, 
                          lapply(1:n_patients, function(i) { 
                            determine_time_steps_since_start(pair$drug_history[i, ])
                          }))


time_steps_since_start

freq_table <- determine_frequency_unique_values(time_steps_since_start, pair$adr_history)

pi0 <- .1
pi1 <- .2
mu <- 5
sigma <- 2

# to make sure that the highest value is indeed 1
normalizing_factor <- dnorm(mu, mu, sigma)


dnorm(time_steps_ago, mu, sigma) / normalizing_factor


param <- c(-1,0,log(5), log(2))

determine_loglikelihood_delayed <- function(param, 
                                               freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  mu_log <- param[3]
  sigma_log <- param[4]
  
  # transform to regular parameters
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  mu <- exp(mu_log)
  sigma <- exp(sigma_log)
  
  # to make sure that the highest value is indeed 1
  normalizing_factor <- dnorm(mu, mu, sigma)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
    #risk_value = exp(-rate * (unique_value - 1))
    risk_value = dnorm(unique_value, mu, sigma) / normalizing_factor
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

determine_loglikelihood_delayed(param, freq_table)
loglikelihood_delayed(param, pair$drug_history, pair$adr_history)

fit <- optim(param,
            expard::loglikelihood_delayed, 
                    freq_table = freq_table,
                    method = "Nelder-Mead",
                    control = list(maxit = 1e8))

beta0 <- fit$par[1]
beta <- fit$par[2]
mu_log <- fit$par[3]
sigma_log <- fit$par[4]

exp(beta0) / (1 + exp(beta0))
exp(beta) / (1 + exp(beta)) 
exp(mu_log)
exp(sigma_log)

freq_table <- expard::determine_frequency_unique_values(time_steps_ago, pair$adr_history)

risk_model <- expard::risk_model_delayed(exp(mu_log), exp(sigma_log))

# extract parameters
beta0 <- param[1]
beta <- param[2]
mu_log <- param[3]
sigma_log <- param[4]

expard::fit_model(pair, model = "delayed")
