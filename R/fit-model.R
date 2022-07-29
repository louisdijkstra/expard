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
#' @section Cohort data oject:
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
#' @param risk_model A risk model function 
#'                    (Default: \code{risk_model_immediate()})
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
fit_model <- function(cohort,
                      risk_model = expard::risk_model_immediate(),
                      start = c(-1,1),
                      method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                                 "Brent"),
                      control = list()) {
  
  # check correctness input 
  cohort <- expard::check_cohort(cohort) 

  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = cohort$n_patients, ncol = cohort$simulation_time) 
  
  # go over all patients 
  for (i in 1:cohort$n_patients) { 
    # go over all timepoints 
    for (t in 1:cohort$simulation_time) { 
      # determine the risk for patient i with drug history 1,2,..,t
      risk <- risk_model(cohort$drug_history[i, 1:t])
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
  class(fit) <- "expardfit3"
  return(fit)
}

#' Print function for the hccd fit_model
#' @export
print.expardfit3 <- function(fit) { 
  cat("expard model fit \n\n")
  cat(sprintf("\t-- no. of patients   : %d\n", fit$n_patients))
  cat(sprintf("\t-- no. of time points: %d\n\n", fit$simulation_time))
  
  if (fit$convergence == 0) {
    cat(green(sprintf("\u2713 CONVERGED\n\n")))
  } else { 
    cat(red(sprintf("\u2717 DID NOT CONVERGE\n\n"))) 
  }
  
  cat(sprintf("Probability of the ADR occurring with\n"))
  cat(sprintf("minimal risk: %.4f\t and \tmaximal risk: %.4f\n", 
              fit$prob_no_adr_with_drug, 
              fit$prob_adr_with_drug))
}
