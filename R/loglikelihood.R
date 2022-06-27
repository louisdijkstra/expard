#' Negative of the Log-Likelihood Function
#' 
#' \code{loglikelihood} returns the negative 
#' log-likelihood function, defined as 
#' \deqn{\pi = \frac{\exp(\beta_0 + \code{risks}\cdot \beta)}{1 + \exp(\beta_0 + \code{risks}\cdot \beta)}} 
#' where \eqn{\pi} is the probability. 
#' 
#' @param betas Vector with two entries: \code{beta0, beta}
#' @param risks Risks determined by the risk model. The risk 
#'              model was applied to the observed drug prescription 
#'              history
#' @export
loglikelihood <- function(betas,
                          risks,
                          adr_history) {

  # extract parameters
  beta0 <- betas[1]
  beta <- betas[2]

  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}
