# This function uses some auxilary functions that can 
# be found at the end of the file


#' Find the Best Model 
#' 
#' Applies all risk models to a given cohort and determines 
#' the AIC, BIC, the test statistic and \eqn{p}-value. 
#' 
#' @param cohort A cohort, see \code{\link{check_cohort}}
#' @param risk_models Risk models that will be applied 
#'                    (Default: \code{c("immediate","immediately_after","extended","long_time_after")})
#' @param params A list with the parameters used for the different risk models
#' @param start,method,control Arguments for the built-in \code{\link{optim}}-function
#' @param verbose Shows progressbar 
#' 
#' @return A \code{tibble} data frame
#' @examples 
#' cohort <- expard::generate_cohort(n_patients = 100, simulation_time = 8)
#' find_best_model(cohort, verbose = FALSE)
find_best_model <- function(cohort,
                            risk_models = c("immediate",
                                            "immediately_after", 
                                            "extended", 
                                            "long_time_after"), 
                            params = list(
                              extension = 1:(cohort$simulation_time-1),
                              duration = 1:(cohort$simulation_time-1), 
                              time_since = 1:(cohort$simulation_time-1)
                            ),
                            start = c(-1,1),
                            method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                                       "Brent"),
                            control = list(),
                            verbose = TRUE) {

  # check correctness input --------------------
  cohort <- hccdanalysis::check_cohort(cohort) 
  check_correctness_input(risk_models, params)
  
  # all results are stored here
  results <- dplyr::tibble()
  
  # set up progress bar ---------------------
  if (verbose) { 
    # determine the number of models that will be fitted 
    n_fits <- determine_number_of_fits(risk_models, params)
    fits_done <- 0 # how many are done
    pb <- txtProgressBar(min = 0, max = n_fits, style = 3)
  }

  
  ### no effect model -----------
  # the no effect model is always fitted. This is needed to determine the 
  # H0 loglikelihood 
  fit <- fit_model(cohort, hccdanalysis::risk_model_no_effect(), start = start, method = method, control = control)
  # update results
  results <- dplyr::bind_rows(results, process_fit_results(fit, "no_effect", 1))
  # set beta and prob_adr_with_drug to NA, since they are not relevant for the 
  # 'no effect' model
  results$beta <- NA
  results$prob_adr_with_drug <- NA
  
  if (verbose) { 
    fits_done <- 1
    setTxtProgressBar(pb, value = fits_done)
  }

  ### immediate model -------------
  if ("immediate" %in% risk_models) { 
    fit <- fit_model(cohort, hccdanalysis::risk_model_immediate(), start = start, method = method, control = control)
    # update results
    results <- dplyr::bind_rows(results, process_fit_results(fit, "immediate", 2))
    
    if (verbose) { 
      fits_done <- fits_done + 1
      setTxtProgressBar(pb, value = fits_done)
    }
  }
  
  ### extended model -------------
  if ("extended" %in% risk_models) { 
    
    # go over all extensions and fit the model for each one 
    for(extension in params$extension) { 
      fit <- fit_model(cohort, hccdanalysis::risk_model_extended(extension), start = start, method = method, control = control)
      # update results
      fit_results <- process_fit_results(fit, "extended", 3)
      fit_results$extension <- extension 
      results <- dplyr::bind_rows(results, fit_results)
      
      if (verbose) { 
        fits_done <- fits_done + 1
        setTxtProgressBar(pb, value = fits_done)
      }
    }
  }
  
  ### immediately after model -----------------------
  if ("immediately_after" %in% risk_models) { 
    
    # go over all extensions and fit the model for each one 
    for(duration in params$duration) { 
      fit <- fit_model(cohort, hccdanalysis::risk_model_immediately_after(duration), start = start, method = method, control = control)
      # update results
      fit_results <- process_fit_results(fit, "immediately_after", 3)
      fit_results$duration <- duration
      results <- dplyr::bind_rows(results, fit_results)
      
      if (verbose) { 
        fits_done <- fits_done + 1
        setTxtProgressBar(pb, value = fits_done)
      }
    }
  }
  
  ### long_time_after model -----------------------
  if ("long_time_after" %in% risk_models) { 
    
    # go over all parameter values
    for (time_since in params$time_since) { 
      for (duration in params$duration) { 
        fit <- fit_model(cohort, hccdanalysis::risk_model_long_time_after(time_since, duration), start = start, method = method, control = control)
        # update results
        fit_results <- process_fit_results(fit, "long_time_after", 4)
        fit_results$duration <- duration
        fit_results$time_since <- time_since
        results <- dplyr::bind_rows(results, fit_results)
      
        if (verbose) { 
          fits_done <- fits_done + 1
          setTxtProgressBar(pb, value = fits_done)
        }
      }
    }
  }
  
  close(pb)

  ### Determine AIC & BIC
  n_observations <- cohort$n_patients * cohort$simulation_time

  # loglikelihood of the no effect model:
  loglikelihood_H0 <- results[results$method == "no_effect",]$loglikelihood

  results <- 
    dplyr::mutate(results,
      AIC = 2*n_params - 2*loglikelihood,
      BIC = log(n_observations)*n_params - 2*loglikelihood,
      test_statistic = -2 * (loglikelihood_H0 - loglikelihood),
      p_value = pchisq(test_statistic, df = n_params - 1, lower.tail = FALSE)
    )
  
  return(results)
}

#' Function for checking the correctness of the input given 
#' to find_best_model
#' 
#' @param risk_models 
#' @param params 
check_correctness_input <- function(risk_models, params) { 
  # check whether the given risk models are correct
  if (length(risk_models) == 0) { 
    stop("no risk model given") 
  }
  
  # check whether the given risk models are valid
  if (!all(risk_models %in% c("immediate","immediately_after", "extended", "long_time_after"))) { 
    stop(sprintf("risk models should be either 'immediate','immediately_after', 'extended' or 'long_time_after'")) 
  }
  
  given_parameters <- names(params)
  if ("extended" %in% risk_models) { 
    if (!("extension" %in% given_parameters)) { 
      stop("The 'extended' model requires the parameter 'extension'") 
    }
  }
  
  if ("immediatedly_after" %in% risk_models) { 
    if (!("duration" %in% given_parameters)) { 
      stop("The 'immediately_after' model requires the parameter 'duration'") 
    }
  }
  
  if ("long_time_after" %in% risk_models) { 
    if ((!("time_since" %in% given_parameters)) || (!("duration" %in% given_parameters))) { 
      stop("The 'long_time_after' model requires the parameters 'time_since' and 'duration'")
    }
  } 
}


#' Auxilary function to process a model fit and return 
#' a tibble 
process_fit_results <- function(fit, method, n_params) { 
  dplyr::tibble(
    method          = method, 
    n_patients      = fit$n_patients, 
    simulation_time = fit$simulation_time, 
    n_params        = n_params, 
    loglikelihood   = fit$loglikelihood, 
    beta0           = fit$beta0, 
    beta            = fit$beta, 
    prob_no_adr_with_drug = fit$prob_no_adr_with_drug, 
    prob_adr_with_drug    = fit$prob_adr_with_drug, 
    convergence     = fit$convergence
  )
}

#' Auxilary function that determines the number of 
#' models that need to be fitted given the input
determine_number_of_fits <- function(risk_models, params) { 
  
  # number of fits. The 'no effect' model is always used
  n_fits <- 1  
  
  # go over the different models and add the number of 
  # model fits needed for each of them
  if ("immediate" %in% risk_models) { n_fits <- n_fits + 1}
  if ("immediately_after" %in% risk_models) { 
    n_fits <- n_fits + length(params$duration) 
  }
  if ("extended" %in% risk_models) { 
    n_fits <- n_fits + length(params$extension) 
  }
  if ("long_time_after" %in% risk_models) { 
    n_fits <- n_fits + length(params$duration) * length(params$time_since) 
  }
  return(n_fits)
}

