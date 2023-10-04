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
  risk_model = c("risk_model_withdrawal(.3)"),
  min_chance_drug = .1,
  avg_duration = 5,
  max_chance_drug = NULL,
  prob_guaranteed_exposed = .1,
  min_chance = .1,
  max_chance = .8,
  verbose = TRUE
)

pair <- cohort[[1]]

# process the data -------------------------------------------------------------

# determine times since last exposure for a drug history of a single patient
determine_time_steps_ago <- function(drug_history) {
  sapply(1:length(drug_history), function(t) { 
    # currently exposed or did not take the drug yet
    if (drug_history[t] == 1 || sum(drug_history[1:t]) == 0) { 
      return(0)  
    } else {
      time_steps_ago = t - max(which(drug_history[1:t] == 1))
      return(time_steps_ago) 
    }
  })
}

n_patients <- nrow(pair$drug_history)

time_steps_ago <- do.call(rbind, 
                          lapply(1:n_patients, function(i) { 
                            determine_time_steps_ago(pair$drug_history[i, ])
                          }))

# determine how often different values from time_steps_ago 
# appear in the data set and how often the ADR occurs for 
# each of them. This is used for speeding up the computation of the 
# loglikeihood
determine_frequency_unique_values <- function(mat, adr_history) {
  unique_values <- sort(unique(as.vector(mat)))
  
  freq_table <-
    lapply(sort(unique_values), function(unique_value) {
      # where does this 'unique_value' occur in the data
      indices <- which(mat == unique_value)
      
      # count how often an ADR occurred at those time points
      n_adr <- sum(adr_history[indices])
      
      c(
        unique_value = unique_value,
        freq       = length(indices),
        # how often does it occur
        n_adr      = n_adr,
        n_no_adr   = length(indices) - n_adr
      )
    })
  
  # create a tibble and return it
  return(as_tibble(do.call(rbind, freq_table)))
}

freq_table <- determine_frequency_unique_values(time_steps_ago, pair$adr_history)
#freq_table <- rename(freq_table, time_since = unique_value)



determine_loglikelihood_withdrawal <- function(param, 
                                               freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  phi <- param[3]
  
  # transform to regular parameters
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  rate <- exp(phi)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
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



param <- c(logit(.1),logit(.2),log(.5))
determine_loglikelihood_withdrawal(param, freq_table)

fit <- optim(c(.1,.1,.1),
             determine_loglikelihood_withdrawal, 
             freq_table = freq_table,
             method = "Nelder-Mead",
             control = list(maxit = 1e8))

r = fit_all_models(pair, models = c("no-association", "current-use", "past-use", "withdrawal"))

beta0 <- fit$par[1]
beta <- fit$par[2]
phi <- fit$par[3]

exp(beta0) / (1 + exp(beta0))
exp(beta) / (1 + exp(beta)) 
exp(phi)


fit2 <- optim(c(.1,.1,.1),
              loglikelihood_withdrawal, 
              freq_table = freq_table,
              method = "Nelder-Mead",
              control = list(maxit = 1e8))

beta0 <- fit2$par[1]
beta <- fit2$par[2]
phi <- fit2$par[3]

exp(beta0) / (1 + exp(beta0))
exp(beta) / (1 + exp(beta)) 
exp(phi)


fit_model(pair, model = "past-use")
fit_model(pair, model = "no-association")

all_fit <- fit_all_models(pair, models = c("no-association", "current-use", "past-use", "withdrawal"))

beta0 <- fit$par[1]
beta <- fit$par[2]
phi <- fit$par[3]

exp(beta0) / (1 + exp(beta0))
exp(beta) / (1 + exp(beta)) 
exp(phi)

determine_loglikelihood_withdrawal(fit$par, freq_table)
expard::loglikelihood_withdrawal(fit$par, pair$drug_history, pair$adr_history)

determine_loglikelihood_withdrawal(c(.1,.1,.1), freq_table)
expard::loglikelihood_withdrawal(c(.1,.1,.1), pair$drug_history, pair$adr_history)


# 'convert' the drug prescriptions. They reflect which period is considered
# to have an increased risk
risks <- matrix(0, nrow = n_patients, ncol = simulation_time)

# go over all patients 
for (i in 1:n_patients) { 
  # go over all timepoints 
  risks[i, ] <- risk_model(pair$drug_history[i, ])
}

tab <- table(risks)

unique_values_risk <- unique(as.vector(risks))

lapply(unique_values_risk, function(value) {
  # find indices in which this values occurs 
  which(risks == value)
})



