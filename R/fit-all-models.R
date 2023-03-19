
# 
# risk_models <- c('no-association', 
#                  'current-use', 
#                  'past-use', 
#                  'withdrawal', 
#                  'delayed',
#                  'decaying', 
#                  'delayed+decaying', 
#                  'long-term')

#' @export
fit_all_models <- function(pair,
                           method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN",
                                      "Brent"),
                           parameters = list()) { 
  
  # models <- c(
  #   'no-association',
  #   'current-use',
  #   'past-use',
  #   'withdrawal',
  #   'delayed',
  #   'decaying',
  #   'delayed+decaying',
  #   'long-term'
  # )
  
  models <- c(
    'no-association',
    'current-use',
    'past-use',
    'withdrawal',
    'delayed',
    'decaying',
    'delayed+decaying',
    'long-term'
  )
  
  res_temp <- lapply(models, function(model) { 
    fit_model(pair, model, method, parameters)  
  })
  
  res <- res_temp[[1]]
  
  for(i in 2:length(res_temp)) { 
    res <- full_join(res, res_temp[[i]])  
  }
  
  do.call("rbind", res_temp)
  
  
  res %>% reduce(full_join)
  
  Reduce(function(...) merge(..., all=T), res)
  
}