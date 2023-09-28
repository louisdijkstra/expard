#' @export
fit_all_models <- function(pair,
                           models = c(
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
  
  res_temp <- lapply(models, function(model) { 
    cat(sprintf("Fitting model %s...\n", model))
    fit_model(pair, model, method, maxiter, parameters)  
  })
  
  res <- res_temp[[1]]
  
  for(i in 2:length(res_temp)) { 
    res <- full_join(res, res_temp[[i]])  
  }
  
  #res$BIC <- res$n_param * log(res$n_patients * res$simulation_time) + 2*res$loglikelihood
  
  # determine the posterior probability for each model 
  # use the parameter setting of the model with the best BIC 
  min_BIC <- min(res$bestBIC)
  
  r <-  res %>% 
    group_by(model) %>% 
    slice_min(n = 1, BIC) %>% 
    mutate(delta_BIC = BIC - min_BIC)
  
  denominator <- sum(exp(-r$delta_BIC/2))
  r <- r %>% mutate(
    posterior = exp(-delta_BIC / 2) / denominator
  ) 
  
  r <- r %>% dplyr::select(model, posterior)
  
  res <- full_join(r, res)
  
  return(res %>% dplyr::arrange(-posterior))
}