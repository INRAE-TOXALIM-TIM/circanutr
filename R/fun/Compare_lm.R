#' Compare Linear Models for Gene Expression
#'
#' @description This function compares different linear models for gene expression data. It fits flat, linear, and circadian models to the data and selects the best model based on the Bayesian Information Criterion (BIC).
#'
#' @param dryList A list containing gene expression results of . It should have a data frame named `results` with gene expression data.
#' @param groupTable A data frame containing sample information, including columns for `sample_name`, `Diet`, and `ZT` (zeitgeber time).
#' @param group A character string specifying the dietary group to subset the data.
#'
#' @return A data frame with BIC values for each model (flat, linear, circadian) and additional columns indicating the best model and comments on the model fit.
#' @export
#'
#' @examples
#' Example usage:
#' Compare_lm(dryList, groupTable, group = "Control")

Compare_lm <- function(dryList, groupTable, group) {
  # Transpose the gene expression data so that genes are columns and samples are rows
  logCPM <- t(dryList$results[, colnames(dryList$results) %in% groupTable$sample_name])
  logCPM <- as.data.frame(logCPM)
  
  # Subset the sample information for the specified dietary group
  pData <- groupTable[, c("sample_name", "Diet", "ZT")]
  
  # Fit flat models (null models with only an intercept)
  flat_models <- apply(logCPM[pData$sample_name[pData$Diet == group], ], 2, function(gene_expr)
    lm(gene_expr ~ 1)) # Model with only an intercept
  
  # Fit linear models
  {
    lin_models <- apply(logCPM[pData$sample_name[pData$Diet == group], ], 2, function(gene_expr)
      lm(gene_expr ~ ZT, data = data.frame(ZT = pData$ZT[pData$Diet == group])))
    
    # Determine the direction of the linear models (up or down)
    models.direction <- sapply(lin_models, function(mdl)
      ifelse(mdl$coefficients[2] < 0, 'down', 'up'))
  }
  
  # Fit circadian models
  circ_models <- apply(logCPM[pData$sample_name[pData$Diet == group], ], 2, function(gene_expr)
    lm(
      gene_expr ~ cos(2 * pi * ZT / 24) + sin(2 * pi * ZT / 24),
      data = data.frame(ZT = pData$ZT[pData$Diet == group])
    ))
  
  # Create a table with BIC values for each model
  {
    BIC_models <- cbind.data.frame(
      flat = sapply(flat_models, BIC),
      linear = sapply(lin_models, BIC),
      circadian = sapply(circ_models, BIC)
    )
    
    # Determine the best model based on the lowest BIC value
    BIC_models$best_model <- as.character(apply(BIC_models, 1, function(x) {
      colnames(BIC_models)[!sapply(x, is.na)][x == min(x, na.rm = T)]
    }))
    
    # Add comments based on the chosen model
    BIC_models$comment <- ifelse(
      dryList$results$chosen_model == 1,
      'flat',
      ifelse(
        dryList$results$chosen_model == 2,
        'loss',
        ifelse(
          dryList$results$chosen_model == 3,
          'gain',
          ifelse(dryList$results$chosen_model == 4, 'sync', 'shift')
        )
      )
    )
    
    # Update comments for linear models with direction information
    BIC_models$comment[BIC_models$best_model == 'linear'] <-
      models.direction[BIC_models$best_model == 'linear']
  }
  
  # Get the R-squared for best model
  BIC_models$Rsqr <- -999
  for (i in 1:length(BIC_models$Rsqr)) {
    BIC_models$Rsqr[i] <- ifelse(
      BIC_models$comment[i] %in% c("flat", "loss"),
      NA,
      ifelse(
        BIC_models$comment[i] %in% c("up", "down"),
        summary(lin_models[[i]])$r.squared,
        ifelse(
          BIC_models$comment[i] %in% c("gain", "sync", "shift"),
          summary(circ_models[[i]])$r.squared, -999)
      )
    )
  }
  
  # Return the data frame with BIC values and model comments
  return(BIC_models)
}


