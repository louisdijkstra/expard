
(drug_history <- generate_drug_history(100, 0.1, .9))

mean(sapply(1:10000, function(i) { 
    sum(generate_drug_history(1000, .00001, 5, guaranteed_exposed = T)) 
  }))

1 / (1 - .1)

(adr_hist <- generate_adr_history(drug_history, expard::risk_model_no_association(), 0, .2))

x <- drug_history
x

past <- 5
simulation_time <- length(x)

sapply(1:simulation_time, function(t) { 
  as.numeric(sum(x[max(1,t-past):t]) != 0)
})

sapply(1:simulation_time, function(t) { 
  as.numeric( any(x[max(1,t-past):t] != 0))
})

get_effect <- function(model, drug_history) { 
  model(drug_history)  
}

get_effect(risk_model_no_association(), drug_history)
get_effect(risk_model_current_use(), drug_history)
get_effect(risk_model_past(3), drug_history)
get_effect(risk_model_past(10), drug_history)
get_effect(risk_model_duration(4), drug_history)

get_effect(risk_model_withdrawal(1), drug_history)

get_effect(risk_model_delayed(3, 1), drug_history)

get_effect(risk_model_decaying(1), drug_history)

plot_risk(risk_model = risk_model_delayed(5, 1))

plot_risk(risk_model = risk_model_decaying(2))

plot_risk(risk_model = risk_model_long_time_after(.01, 300))

plot_risk(drug_history = c(rep(0,4), rep(1,10), rep(0,10)), risk_model = risk_model_delayed_decaying(5, 2, .3))


