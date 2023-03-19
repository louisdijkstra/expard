#' Estimate
#' 
#' \code{estimate} fits a specified risk model \code{risk_model}. 
#' The maximum likelihood is determined numerically. 
#' 
#' @section Pair data oject:
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
#' @param pair The data for a single drug-ADR pair. See the \code{\link{generate_cohort}} 
#'             function. A pair is a list of two matrices, see below
#' @param risk_model A risk model function 
#'                    (Default: \code{risk_model_current_use()})
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
#' @seealso \code{\link{generate_cohort}}
#' @examples 
#' cohort <- expard::generate_cohort()
#' 
#' fit <- expard::estimate(cohort[[1]], risk_model = expard::risk_model_current_use())
#' @export
estimate <- function(pair,
                     risk_model = expard::risk_model_current_use(),
                     start = c(-1, 1),
                     method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                                "Brent"),
                     control = list()) {
  
  
  
  
  n_patients <- nrow(pair$drug_history)
  simulation_time <- ncol(pair$drug_history)
  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
  
  # go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(pair$drug_history[i, ])
  }

  # find the optimal values for beta0 and beta
  res <- optim(start,
               expard::loglikelihood, 
               risks = risks,
               adr_history = pair$adr_history,
               method = method[1],
               control = list())

  beta0 <- res$par[1]
  beta <- res$par[2]
  
  fit <- list(
    n_patients = n_patients, 
    simulation_time = simulation_time,
    loglikelihood = res$value,
    beta0 = beta0,
    beta = beta,
    prob_no_adr_with_drug = exp(beta0) / (1 + exp(beta0)),
    prob_adr_with_drug = exp(beta0 + beta) / (1 + exp(beta0 + beta)),
    converged = res$convergence
  )
  class(fit) <- "expardfit"
  return(fit)
}

#' Print function for the expardfit
#' @export
print.expardfit <- function(fit) { 
  cat("expard model fit \n\n")
  cat(sprintf("\t-- no. of patients   : %d\n", fit$n_patients))
  cat(sprintf("\t-- no. of time points: %d\n\n", fit$simulation_time))
  
  if (fit$converged == 0) {
    cat(green(sprintf("\u2713 CONVERGED\n\n")))
  } else { 
    cat(red(sprintf("\u2717 DID NOT CONVERGE\n\n"))) 
  }
  
  cat(sprintf("Log-likelihood: %.4f\n\n", fit$loglikelihood))
  cat(sprintf("Probability of the ADR occurring with\n"))
  cat(sprintf("minimal risk: %.4f\t and \tmaximal risk: %.4f\n", 
              fit$prob_no_adr_with_drug, 
              fit$prob_adr_with_drug))
}
