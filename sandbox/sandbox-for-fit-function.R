library(expard)

cohort <- generate_cohort(n_patients = 20, risk_model = rep("risk_model_delayed(3, 1)", 4), simulation_time = 10, n_drug_ADR_pairs = 4, )
tables <- create2x2tables(cohort)

expard::create2x2tables(cohort)

m <- fit_model2(cohort[[1]], model = 'past-use')
m0 <- fit_model2(cohort[[1]], model = "current-use")

f <- estimate(cohort[[1]], risk_model = expard::risk_model_decaying(4))

res <- data.frame(mu = seq(1, 5, by = 0.1))

res$logl <- sapply(res$mu, function(rate) { 
  f <- fit_model(cohort[[1]], risk_model = expard::risk_model_delayed(rate, 1)) 
  f$loglikelihood
})

plot(res)


table <- expard::create2x2table(cohort[[1]], method = "patient")
tables <- expard::create2x2tables(cohort)

f2 <- fit_model2(cohort[[1]], model = "no-association")
fit_model2(cohort[[1]], model = "current-use")


f <- fit_model(cohort[[1]], risk_model = expard::risk_model_past(4))

f2 <- fit_model2(cohort[[1]], model = "past-use", parameters = list(past = 4))


pair <- generate_drug_ADR_pair()




cohort <-
  generate_cohort(
    n_patients = 1000,
    risk_model = "risk_model_current_use()",
    simulation_time = 20,
    n_drug_ADR_pairs = 1,
    min_chance = .1, 
    max_chance = .8
  )

fit_model2(cohort[[1]], model = "current-use")
