estimates_table <-  function(model, row_names) {
  
  # get model summary
  model_summary <- summary(model)
  
  # get number of coefficients
  num_coeff <-  nrow(model_summary$coefficients)
  
  # extract estimates and confidence intervals for each predictor
  estimates <- model_summary$coefficients[c(1:num_coeff), c("Estimate", "Std. Error")]
  colnames(estimates) <- c("beta", "SE")
  estimates <- as.data.frame.array(estimates)
  
  # calculate confidence intervals and t-values
  estimates$CI_low <- round(estimates$beta - 1.96 * estimates$SE, digits = 2)
  estimates$CI_high <- round(estimates$beta + 1.96 * estimates$SE, digits = 2)
  estimates$t_value <- estimates$beta / estimates$SE
  estimates$p_value <- round(2 * pt(abs(estimates$t_value), df.residual(model),
                                    lower.tail = FALSE), digits = 3)
  estimates$SE <- NULL
  estimates$beta <- round(estimates$beta, digits = 2)
  estimates$t_value  <- round(estimates$t_value, digits = 2)
  
  # set column names
  colnames(estimates) <- c("$\\beta$", "CI low", "CI high", "t-value", "p-value")
  rownames(estimates) <- row_names
  
  # create LaTeX table
  table_latex <- xtable(estimates, align = c("l", "r", "r", "r", "r", "r"))
  # table_latex <- print.xtable(table_latex)
  
  # cat(table_latex)
  print(estimates)
  
}