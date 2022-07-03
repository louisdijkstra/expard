library(expard)

risk_model_temp <- function() {
  fn <- function(drug_history, ...) { 
    if (drug_history[length(drug_history)]) { 
      1  # highest risk
    } else {
      0  # lowest risk
    }
  }

  model <- list(
    name = "immediate",
    get_risk = fn,
    n_params = 2,
    params = c("beta0", "beta"),
    constraints = c(NA, NA)
  )
  
  class(model) <- "expard-model"
  return(model)
}

# type of constraints: 
# discrete or continuous
# NA, c(NA, 4), c(0, NA), c(0, 1)

r = risk_model_temp()
r$fn(c(0,0,1,1))


