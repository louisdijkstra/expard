library(expard)

cohort <- generate_cohort(n_patients = 5, risk_model = rep("risk_model_delayed(3, 1)", 4), simulation_time = 5, n_drug_ADR_pairs = 4, )
tables <- create2x2tables(cohort)

expard::create2x2tables(cohort)

m <- fit_model(cohort[[1]], model = 'past-use')
m <- fit_model(cohort[[1]], model = 'no-association')
m <- fit_model(cohort[[1]], model = 'withdrawal')
m <- fit_model(cohort[[1]], model = 'delayed')
m <- fit_model(cohort[[1]], model = 'decaying')
m <- fit_model(cohort[[1]], model = 'delayed+decaying')
m <- fit_model(cohort[[1]], model = 'long-term')


r = fit_all_models(cohort[[1]])

model = c(
  'no-association',
  'current-use',
  'past-use',
  'withdrawal',
  'delayed',
  'decaying',
  'delayed+decaying',
  'long-term'
)

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

param = c(.1, .2, 2)

drug_history = pair$drug_history
adr_history = pair$adr_history




pair <- cohort[[1]]

param = c(.1, .2, 2)
loglikelihood_withdrawal(param, pair$drug_history, pair$adr_history)

param = c(.1, .8, 0, 0)
loglikelihood_delayed(param, pair$drug_history, pair$adr_history)


# find the optimal values for beta0 and beta
res <- optim(c(0,0,-1),
             expard::loglikelihood_withdrawal, 
             drug_history = pair$drug_history,
             adr_history = pair$adr_history,
             method = "L-BFGS",
             control = list())


lower=c(0, -Inf, -Inf, 0), upper=rep(Inf, 4),
method="L-BFGS-B"

fit_model2(pair, model = "current-use")
fit_model2(pair, model = "withdrawal")



cohort <-
  generate_cohort(
    n_patients = 100,
    risk_model = "risk_model_decaying(.5)",
    simulation_time = 20,
    n_drug_ADR_pairs = 1,
    min_chance = .1, 
    max_chance = .2
  )
pair <- cohort[[1]]
#fit_model(pair, model = "no-association")

res = fit_all_models(pair)



model = c(
  'no-association',
  'current-use',
  'past-use',
  'withdrawal',
  'delayed',
  'decaying',
  'delayed+decaying',
  'long-term'
)
