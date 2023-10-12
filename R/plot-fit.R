#' Plot the Fit 
#' 
#' Plots the BIC for the different risk models that were fitted by the 
#' function \code{\link{fit_all_models}}. 
#' 
#' @param fit The fit 
#' @param xlab The x-label
#' @param past If not \code{NULL}, a select number of values for past are 
#'             used for the plot
#' 
#' @return ggplot 
#' @export
plot_fit <- function(fit, 
                     xlab = "model", 
                     past_values = NULL) {
  
  
  
  if('past-use(1)' %in% fit$model) {
    simulation_time <- fit$simulation_time[1]
    
    past_models <- sprintf('past-use(%d)', 1:(simulation_time-1))
    
    fit$model <- sapply(1:nrow(fit), function(i) {
      if(grepl('past-use', fit$model[i], fixed = TRUE)) {
        return('past-use')
      } else {
        fit$model[i]
      }
    })
  }
  
  if (!is.null(past_values)) {
    fit <- fit %>% filter(model != 'past-use' | past %in% past_values)
  }

  best_fit <- fit %>% group_by(model) %>% 
    filter(BIC == min(BIC))  %>% 
    arrange(BIC)
  
  best_fit <- fit %>% group_by(model) %>% 
    filter(BIC == min(BIC))  %>% 
    arrange(BIC)
  
  min_BIC <- min(best_fit$BIC)
  max_BIC <- max(best_fit$BIC)
  
  ggplot(best_fit) +
    geom_bar(aes(x = reorder(model, BIC), y = BIC), stat="identity") + 
    coord_cartesian(ylim=c(min_BIC,max_BIC)) + 
    scale_y_continuous(expand = expansion(mult = c(0.1, .1))) + #expand = c(0, 100)) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("model") 
}
