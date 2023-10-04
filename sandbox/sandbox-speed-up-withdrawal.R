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
  risk_model = c("risk_model_withdrawal(.5)"),
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

# determine the time_steps_ago for every patient
# time_steps_ago <- unlist(
#   lapply(1:n_patients, function(i) { 
#    determine_time_steps_ago(pair$drug_history[i, ])
#     }))

time_steps_ago <- do.call(rbind, 
  lapply(1:n_patients, function(i) { 
    determine_time_steps_ago(pair$drug_history[i, ])
    }))

# dim(time_steps_ago)

# determine the frequencies with which time_steps_ago appear in the data set
# and how often the ADR occurs etc.

# unique(as.vector(time_steps_ago))

freq_table <- lapply(sort(unique(as.vector(time_steps_ago))), function(time_since) {
  
  # where does this 'time_since' occur in the data
  indices <- which(time_steps_ago == time_since)
  
  # count how often an ADR occurred at those time points
  n_adr <- sum(pair$adr_history[indices])
  
  c(
    time_since = time_since, 
    freq       = length(indices),  # how often does it occur
    n_adr      = n_adr, 
    n_no_adr   = length(indices) - n_adr
  )
})

# turn into a tibble
freq_table <- as_tibble(do.call(rbind, freq_table))


# freq_table$freq %>% sum()
# freq_table$n_adr %>% sum()
# sum(pair$adr_history)

determine_loglikelihood_withdrawal <- function(param, 
                                               freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  phi <- param[3]
  
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  rate <- exp(phi)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
    risk_value = exp(-rate * (time_since - 1))
  )
  
  # the risk when the time_since is 0 (never exposed or currently exposed)
  # is zero
  freq_table$risk_value[1] <- 0
  
  # p = exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  
  # determine the log-likelihood given the risks
  freq_table <- freq_table %>% rowwise() %>% 
    mutate(
      #P = exp(beta0 + risk_value * beta) / (1 + exp(beta0 + risk_value * beta)),
      #loglikelihood = -1 * n_adr * log(P) + n_no_adr*log(1 - P))
      #
      loglikelihood = -1*(n_adr*log((pi1 - pi0)*risk_value + pi0) + 
                                       n_no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))
      #loglikelihood = -1*(n_adr*log((pi1 - pi0)*risk_value + pi0) + 
      #                                                     n_no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))
  #sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
  sum(freq_table$loglikelihood)
}



param <- c(logit(.1),logit(.2),log(.5))
determine_loglikelihood_withdrawal(param, freq_table)

fit <- optim(c(.1,.1,.1),
             determine_loglikelihood_withdrawal, 
             freq_table = freq_table,
             method = "Nelder-Mead",
             control = list(maxit = 1e8))

beta0 <- fit$par[1]
beta <- fit$par[2]
phi <- fit$par[3]

exp(beta0) / (1 + exp(beta0))
exp(beta0 + beta) / (1 + exp(beta0 + beta))

exp(beta0) / (1 + exp(beta0))
exp(beta) / (1 + exp(beta)) 
exp(phi)


fit2 <- optim(c(.1,.1,.1),
             loglikelihood_withdrawal, 
             drug_history = pair$drug_history, 
             adr_history = pair$adr_history,
             method = "Nelder-Mead",
             control = list(maxit = 1e8))

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







simulation_time <- ncol(pair$drug_history)

# determine times since last exposure

determine_time_steps_ago <- function(drug_history) {
  sapply(1:simulation_time, function(t) { 
    # currently exposed or did not take the drug yet
    if (drug_history[t] == 1 || sum(drug_history[1:t]) == 0) { 
      return(0)  
    } else {
      time_steps_ago = t - max(which(drug_history[1:t] == 1))
      return(time_steps_ago) 
    }
  })
}


time_steps_ago <- unlist(lapply(1:n_patients, function(i) { 
  determine_time_steps_ago(pair$drug_history[i, ])
  })
)

freq <- table(time_steps_ago)

freq <- list(
  time_since = as.integer(names(freq)), 
  frequency = freq
)

res <- lapply(freq$time_since, function(time_since) {
  indices <- which(time_steps_ago == time_since)
  n_adr = sum(pair$adr_history[indices])
  c(
    time_since = time_since, 
    freq = length(indices),
    n_adr  = n_adr, 
    no_adr = length(indices) - n_adr
  )
})

res <- as_tibble(do.call(rbind, res))

res$risk_value[1] <- 0





freq_table <-  res

res$risk_value <- exp(-rate * (unlist(res$time_since) - 1))

pi1 <- .1 
pi0 <- .05

res <- res %>% mutate(loglikelihood = -1*(n_adr * log((pi1 - pi0)*risk_value + pi0) + 
                 no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))

sum(res$loglikelihood)

#fit$p0 = exp(beta0) / (1 + exp(beta0))
#fit$p1 = exp(beta0 + beta) / (1 + exp(beta0 + beta))
#fit$rate <- exp(res$par[3])

logit <- function(p) {
  log(p / (1 - p))
}

loglikelihood_withdrawal(c(logit(.05),logit(.4),log(rate)), pair$drug_history, pair$adr_history)


my_loglikelihood_withdrawal <- function(param, 
                                        freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  phi <- param[3]
  
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  rate <- exp(phi)

  freq_table <- freq_table %>% mutate(
    risk_value = exp(-rate * (time_since - 1))
  )
  
  freq_table$risk_value[1] <- 0
  
  freq_table <- freq_table %>% mutate(loglikelihood = -1*(n_adr * log((pi1 - pi0)*risk_value + pi0) + 
                                              no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))
  
  sum(freq_table$loglikelihood)
}

param <- c(logit(.1), logit(.2), log(2))

maxiter = 10000000
fit <- optim(c(logit(.01),logit(.02),log(1)),
             my_loglikelihood_withdrawal, 
             freq_table = freq_table,
             method = "Nelder-Mead",
             control = list(maxit = maxiter))

beta0 <- fit$par[1]
beta <- fit$par[2]

exp(fit$par[1]) / (1 + exp(fit$par[1]))
exp(fit$par[2]) / (1 + exp(fit$par[2]))
exp(fit$par[3])


 exp(beta0) / (1 + exp(beta0))
exp(beta0 + beta) / (1 + exp(beta0 + beta))
log(fit$par[3])

exp(beta0) / (1 + exp(beta0))
exp(beta) / (1 + exp(beta))
exp(fit$par[3])

fit_model(pair, model = "withdrawal")
# 
# fit$p1 = exp(beta0 + beta) / (1 + exp(beta0 + beta))
# fit$rate <- exp(res$par[3])
# 
# my_loglikelihood_withdrawal(c(logit(pi0), logit(pi1), log(rate)), res)
# 
# 
# 
# loglikelihood_withdrawal <- function(param,
#                                      drug_history,
#                                      adr_history) {
#   
#   # extract parameters
#   beta0 <- param[1]
#   beta <- param[2]
#   phi <- param[3]
#   
#   risk_model <- expard::risk_model_withdrawal(exp(phi))
#   
#   
#   n_patients <- nrow(drug_history)
#   simulation_time <- ncol(drug_history)
#   
#   # 'convert' the drug prescriptions. They reflect which period is considered
#   # to have an increased risk
#   risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
#   
#   # go over all patients 
#   for (i in 1:n_patients) { 
#     # go over all timepoints 
#     risks[i, ] <- risk_model(drug_history[i, ])
#   }
#   
#   #risks <- Matrix(apply(drug_history, 1, function(d) risk_model(d)))
#   
#   P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
#   -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
# }
# 
# 
# logl <- n_adr * log((pi1 - pi0)*risk_value + pi0) + 
#   no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)
# 
# 
# 
# freq <- list(time_since = as.integer(names(x)), 
#              freq = x)
# 
# rate = 5
# determine_loglikelihood_withdrawal <- function(rate, time_since, freq) { 
#     risks <- exp(-rate * (time_since - 1))
# }
# 
# 
# exp(-rate * (freq$time_since - 1))
# 
# 
# 
# lapply(1:n_patients, function(i) { 
#   # go over all timepoints 
#   timepoints_since <- rep(0,simulation_time)
#   drug_history <- pair$drug_history[i, ]
#   adr_history <- pair$adr_history[i, ]
#   sapply(2:simulation_time, function(t) {
#     if (drug_history[t] == 0) { 
#       if (timepoints_since[t-1] != 0) {
#         time_points_since[t] <<- timepoints_since[t-1] + 1 
#       } else {
#         timepoints_since[t] <<- 0
#       }
#     }
#   })
#   
#   timepoints_since
# })
# 
# 
# 
# 
# risks[9943]
# 
# 
