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
#' @param maxiter Maximum number iterations for the \code{\link{optim}}-solver
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
fit_model <- function(pair,
                      model = c(
                        'no-association',
                        'current-use',
                        'past-use',
                        'withdrawal',
                        'delayed',
                        'decaying',
                        'delayed+decaying',
                        'long-term'
                      ),
                      method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN",
                                 "Brent"),
                      maxiter = 1000, 
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
    
    fit$p <- pi
    fit$n_param <- 1
    fit$loglikelihood <- -1*((table$a + table$b) * log(pi) + (table$c + table$d) * log(1 - pi))
    fit$converged <- TRUE
    fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
    fit$bestBIC <- min(fit$BIC)
    
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
    
    fit$p1 <- pi1
    fit$p0 <- pi0
    
    fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
    fit$bestBIC <- min(fit$BIC)
    
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
  
  
  
  if (model[1] == "withdrawal") { 

    res <- optim(c(0,0,-1),
                 expard::loglikelihood_withdrawal, 
                 drug_history = pair$drug_history,
                 adr_history = pair$adr_history,
                 method = "Nelder-Mead",
                 control = list(maxit = maxiter))
    
    beta0 <- res$par[1]
    beta <- res$par[2]
    fit$p0 = exp(beta0) / (1 + exp(beta0))
    fit$p1 = exp(beta0 + beta) / (1 + exp(beta0 + beta))
    fit$rate <- exp(res$par[3])
    
    fit$n_param <- 3
    fit$loglikelihood <-res$value
    fit$converged <- res$convergence == 0
    
    fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
    fit$bestBIC <- min(fit$BIC)
    
    return(fit)
  }
  
  
  if (model[1] == "delayed") { 
    
    res <- optim(c(0,0,0,0),
                 expard::loglikelihood_delayed, 
                 drug_history = pair$drug_history,
                 adr_history = pair$adr_history,
                 method = "Nelder-Mead",
                 control = list(maxit = maxiter))
    
    beta0 <- res$par[1]
    beta <- res$par[2]
    fit$p0 = exp(beta0) / (1 + exp(beta0))
    fit$p1 = exp(beta0 + beta) / (1 + exp(beta0 + beta))
    fit$mu <- exp(res$par[3])
    fit$sigma <- exp(res$par[4])
    
    fit$n_param <- 4
    fit$loglikelihood <- res$value
    fit$converged <- res$convergence == 0
    
    fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
    fit$bestBIC <- min(fit$BIC)
    
    return(fit)
  }
  

  if (model[1] == "decaying") { 
    
    res <- optim(c(0,0,-1),
                 expard::loglikelihood_decaying, 
                 drug_history = pair$drug_history,
                 adr_history = pair$adr_history,
                 method = "Nelder-Mead",
                 control = list(maxit = maxiter))
    
    beta0 <- res$par[1]
    beta <- res$par[2]
    fit$p0 = exp(beta0) / (1 + exp(beta0))
    fit$p1 = exp(beta0 + beta) / (1 + exp(beta0 + beta))
    fit$rate <- exp(res$par[3])
    
    fit$n_param <- 3
    fit$loglikelihood <- res$value
    fit$converged <- res$convergence == 0
    
    fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
    fit$bestBIC <- min(fit$BIC)
    
    return(fit)
  }
  
  
  
  if (model[1] == "delayed+decaying") { 
    
    res <- optim(c(0,0,0,0,-1),
                 expard::loglikelihood_delayed_decaying, 
                 drug_history = pair$drug_history,
                 adr_history = pair$adr_history,
                 method = "Nelder-Mead",
                 control = list(maxit = maxiter))
    
    beta0 <- res$par[1]
    beta <- res$par[2]
    fit$p0 = exp(beta0) / (1 + exp(beta0))
    fit$p1 = exp(beta0 + beta) / (1 + exp(beta0 + beta))
    fit$mu <- exp(res$par[3])
    fit$sigma <- exp(res$par[4])
    fit$rate <- exp(res$par[5])
    
    fit$n_param <- 5
    fit$loglikelihood <- res$value
    fit$converged <- res$convergence == 0
    
    fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
    fit$bestBIC <- min(fit$BIC)
    
    return(fit)
  }
  
  
  
  if (model[1] == "long-term") { 
    
    res <- optim(c(0,0,0,0),
                 expard::loglikelihood_long_term, 
                 drug_history = pair$drug_history,
                 adr_history = pair$adr_history,
                 method = "Nelder-Mead",
                 control = list(maxit = maxiter))
    
    beta0 <- res$par[1]
    beta <- res$par[2]
    fit$rate <- exp(res$par[3])
    fit$delay <- exp(res$par[4])
    
    fit$n_param <- 4
    fit$loglikelihood <- res$value
    fit$converged <- res$convergence == 0
    
    fit$BIC <- fit$n_param * log(fit$n_patients * fit$simulation_time) + 2 * fit$loglikelihood
    fit$bestBIC <- min(fit$BIC)
    
    return(fit)
  }
  
  stop("model not known")
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
