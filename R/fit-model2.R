#' Fit a Model to Cohort Data
#' 
#' \code{fit_model} fits a specified risk model \code{risk_model}
#' to cohort data. 
#' 
#' @section Patient models: 
#' \code{fit_model} can currently only deal with the 
#' \code{\link{patient_model_uninformative}}. More complex
#' patient models, such as \code{\link{patient_model_sex}} 
#' require the estimation of more parameters which is currently
#' not supported. 
#' @section Cohort data object:
#' Minimal requirement is that the cohort is a list with 
#' two matrix with the items
#' \itemize{
#'   \item{\code{drug_history}  A binary matrix of size 
#'        \code{n_patients x simulation_time}  describing
#'        the drug prescriptions.}
#'   \item{\code{adr_history}  A binary matrix of size 
#'        \code{n_patients x simulation_time} describing
#'        ADR histories}
#' }
#' See \code{\link{check_cohort}} and \code{\link{generate_cohort}}
#' for more details on the \code{cohort} object. 
#' 
#' @param cohort A cohort dataset. See details below. 
#' @param model Label for a risk model. Can be either ...
#' @param start Starting point for the base \code{\link{optim}}-solver
#'              (Default: \code{c(-1,1)})
#' @param method Methods used by the base \code{\link{optim}}-solver
#' @param control List used by the base \code{\link{optim}}-solver
#'               
#' @return A model fit. A list with the items 
#'       \item{\code{n_patients}}{Total number of patients in the dataset}
#'       \item{\code{simulation_time}}{Total number of time points}
#'       \item{\code{loglikelihood}}{The maximal log-likelihood value}
#'       \item{\code{beta0}}{Intercept for the logistic model}
#'       \item{\code{beta}}{Parameter for the logistic model}
#'       \item{\code{prob_no_adr_with_drug}}{Estimated probability
#'                          of no ADR occurring}
#'       \item{\code{prob_ADR_with_drug}}{Estimated probability
#'                          of ADR occurring}
#'       \item{\code{convergence}}{Convergence value of \code{\link{optim}}. 
#'                \code{0} means that the algorithm converged}
#'                          
#' @seealso \code{\link{check_cohort}},\code{\link{generate_cohort}}
#' @examples 
#' cohort <- expard::generate_cohort(n_patients = 1000, 
#'                                   simulation_time = 100, 
#'                                   risk_model = expard::risk_model_immediate(), 
#'                                   verbose = TRUE, 
#'                                   min_chance_drug = probability_model_constant(.3), 
#'                                   max_chance_drug = probability_model_constant(.7), 
#'                                   min_chance_adr = probability_model_constant(.3),
#'                                   max_chance_adr = probability_model_constant(.6))
#' # fit the no effect model
#' fit_model(cohort, risk_model = expard::risk_model_no_effect()) 
#' 
#' # fit the true immediate effect model
#' fit_model(cohort, risk_model = expard::risk_model_immediate())
#' # note that the estimators are close to the truth (.3 and .6) 
#' @export
fit_model2 <- function(pair,
                       model = c('no-association', 
                                 'current-use', 
                                 'past-use', 
                                 'withdrawal', 
                                 'delayed',
                                 'decaying', 
                                 'delayed+decaying', 
                                 'long-term'),
                       method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN",
                                 "Brent"),
                       parameters = list()) {
 
  # initialize the fit --------------------------------
  fit <- data.frame(
    n_patients = nrow(pair$drug_history), 
    simulation_time = ncol(pair$drug_history), 
    model = model[1], 
    n_param = NA, 
    loglikelihood = NA, 
    converged = NA
  )
  
  #class(fit) <- c(class(fit), "expardmodel")
  
  # No association model -------------------------------------------------------
  if (model[1] == "no-association") { 
    
    # determine the 2x2 table
    table <- expard::create2x2table(pair, method = "time-point")
    
    pi <- (table$a + table$b) / table$n
    
    fit$pi <- pi
    fit$n_param <- 1
    fit$loglikelihood <- -1*((table$a + table$b) * log(pi) + (table$c + table$d) * log(1 - pi))
    fit$converged <- TRUE
    
    return(fit)
  }
  
  
  if (model[1] == "current-use") {
    # create 2x2 tables 
    table <- expard::create2x2table(pair, method = "time-point")
    
    pi1 <- table$a / (table$a + table$c) 
    pi0 <- table$b / (table$b + table$d)
    
    fit$n_param <- 2
    fit$loglikelihood <- -1*(table$a)*log(pi1) - (table$c)*log(1 - pi1) - table$b*log(pi0) - (table$d)*log(1 - pi0)
    fit$converged <- TRUE
    
    fit$pi1 <- pi1
    fit$pi0 <- pi0
    
    return(fit)
  }
  
  
  
  if (model[1] == "past-use") { 
    
    past <- 1:(fit$simulation_time - 1)
    
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
    
    estimates <- lapply(past, function(d) { 
      res <- expard::estimate(pair, 
                       risk_model = expard::risk_model_past(d))
      res$past <- d
      return(res)
      })
    
    fit$loglikelihood <- sapply(estimates, function(est) est$loglikelihood)
    fit$pi0 <- sapply(estimates, function(est) est$adr_while_not_at_risk / (est$adr_while_not_at_risk + est$no_adr_while_not_at_risk))
    fit$pi1 <- sapply(estimates, function(est) est$adr_while_at_risk / (est$adr_while_at_risk + est$no_adr_while_at_risk))
    fit$converged <- TRUE

    # select the best model
    #best <- estimates[[1]]
    
    #lapply(estimates, function(est) { 
    #    if (est$loglikelihood < best$loglikelihood) { 
    #      best <<- est   
    #    }
    #  })
    
    # fit$est <- list(
    #   pi0 = best$prob_no_adr_with_drug, 
    #   pi1 = best$prob_adr_with_drug, 
    #   past = best$past
    # )
    # 
    # fit$n_param <- 3
    # fit$loglikelihood <- best$loglikelihood
    # fit$converged <- best$converged
    
    return(fit)
  }
  
  
  if (model[1] == "'withdrawal'") { 
    
    # number of time points in the past
    

    pi0 <- adr_while_not_at_risk / (adr_while_not_at_risk + no_adr_while_not_at_risk)
    pi1 <- adr_while_at_risk / (adr_while_at_risk + no_adr_while_at_risk)
    
    fit$est <- list(
      pi1 = pi1, 
      pi0 = pi0
    )
    
    fit$n_param <- 3
    fit$loglikelihood <- -1*adr_while_at_risk*log(est$pi1) - 
                      no_adr_while_at_risk*log(1 - est$pi1) - 
                      adr_while_not_at_risk*log(est$pi0) - 
                      no_adr_while_not_at_risk*log(1 - est$pi0)
    fit$converged <- TRUE
    
    return(fit)
  }

  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = cohort$n_patients, ncol = cohort$simulation_time) 
  
  # go over all patients 
  for (i in 1:cohort$n_patients) { 
    # go over all timepoints 
    for (t in 1:cohort$simulation_time) { 
      # determine the risk for patient i with drug history 1,2,..,t
      risk <- risk_model$fn(cohort$drug_history[i, 1:t])
      risks[i,t] <- risk
    }
  }

  # find the optimal values for beta0 and beta
  res <- optim(start,
               loglikelihood, 
               risks = risks,
               adr_history = cohort$adr_history,
               method = method[1],
               control = list())

  beta0 <- res$par[1]
  beta <- res$par[2]

  fit <- list(
    n_patients = cohort$n_patients, 
    simulation_time = cohort$simulation_time,
    loglikelihood = -1 * res$value,
    beta0 = beta0,
    beta = beta,
    prob_no_adr_with_drug = exp(beta0) / (1 + exp(beta0)),
    prob_adr_with_drug = exp(beta0 + beta) / (1 + exp(beta0 + beta)),
    convergence = res$convergence
  )
  class(fit) <- "expardfit"
  return(fit)
}

#' Print function for the hccd fit_model
#' @export
print.expardfit2 <- function(fit) { 
  cat("expard model fit\n")
  cat("----------------\n\n")
  cat(sprintf("no. of patients    : %d\n", fit$n_patients))
  cat(sprintf("no. of time points : %d\n\n", fit$simulation_time))
  
  if (fit$converged) {
    cat(green(sprintf("\u2713 CONVERGED\n\n")))
  } else { 
    cat(red(sprintf("\u2717 DID NOT CONVERGE\n\n"))) 
  }
  
  cat(sprintf("negative log-likelihood : %g\n", fit$loglikelihood))
  cat(sprintf("no. of parameters       : %d\n\n", fit$n_param))
  
  cat("estimates:\n\n")
  for (i in 1:fit$n_param) { 
     cat(sprintf("\t-- %s : %g\n", names(fit$est)[i], fit$est[i]))
  }
}
