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
fit_model2 <- function(cohort,
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
 
  
  
  
  if (model[1] == "no_association") { 
    # create 2x2 tables 
    tables <- expard::create2x2tables(cohort, method = "time-point")
    
    est <- lapply(tables, function(table) { 
      list(pi = (table$a + table$b) / table$n)
    })
    
    logl <- lapply(1:length(tables), function(i) {
      table <- tables[[i]]
      -1*(table$a + table$b) * log(est[[i]]$pi) - (table$c + table$d) * log(1 - est[[i]]$pi)
    })
      
    converged <- rep(TRUE, length(tables))
  }
  
  if (model[1] == "current_use") {
    # create 2x2 tables 
    table <- expard::create2x2table(cohort, method = "time-point")
    
    est <- list(
        pi0 = table$b / (table$b + table$d),
        pi1 = table$a / (table$a + table$c)
      ) 
    loglikelihood <- -1*(table$a)*log(est$pi1) - (table$c)*log(1 - est$pi1) - table$b*log(est$pi0) - (table$d)*log(1 - est$pi0)
    converged <- TRUE
  }
  
  if (model[1] == "past") { 
    
    # number of time points in the past
    k <- parameters$k 
    
    if (k >= cohort$simulation_time) { 
      stop(sprintf("past risk model - number of time points in the past k = %d is larger than the simulation time T = %d", k, simulation_time)) 
    }
    
    # number of times the ADR does or does not happen
    # while being at full risk or not
    adr_while_at_risk        <- 0 
    no_adr_while_at_risk     <- 0 
    adr_while_not_at_risk    <- 0 
    no_adr_while_not_at_risk <- 0 
    
    for(i in 1:cohort$n_patients) { 
      # go over all time-points
      for (t in (k+1):cohort$simulation_time) { 
        
        # if ADR happens at time point t for patient k
        if (cohort$adr_history[i,t]) { 
          if (any(cohort$drug_history[i, 1:t] == 1)) # took drug in the last k time points
            adr_while_at_risk <- adr_while_at_risk + 1
          else { 
            adr_while_not_at_risk <- adr_while_not_at_risk + 1
          }
        } else { # ADR did not occur
          if (any(cohort$drug_history[i, 1:t] == 1)) # took drug in the last k time points
            no_adr_while_at_risk <- no_adr_while_at_risk + 1
          else { 
            no_adr_while_not_at_risk <- no_adr_while_not_at_risk + 1
          }
        }
        
      }
    }

    est <- list(
      pi0 = adr_while_not_at_risk / (adr_while_not_at_risk + no_adr_while_not_at_risk),
      pi1 = adr_while_at_risk / (adr_while_at_risk + no_adr_while_at_risk)
    ) 
    
    loglikelihood <- -1*adr_while_at_risk*log(est$pi1) - 
                      no_adr_while_at_risk*log(1 - est$pi1) - 
                      adr_while_not_at_risk*log(est$pi0) - 
                      no_adr_while_not_at_risk*log(1 - est$pi0)
    converged <- TRUE
  }

  fit <- list(
    est = est, 
    n_param = length(est[[1]]),
    loglikelihood = logl, 
    n_patients = cohort$n_patients, 
    simulation_time = cohort$simulation_time,
    converged = converged,
    parameters = parameters
  )
  class(fit) <- "expardfit"
  return(fit)
  
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
print.expardfit <- function(fit) { 
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
