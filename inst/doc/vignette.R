## ----setup, include = FALSE, fig.align='center', warning = F, message=F-------
options(tinytex.verbose = TRUE)
knitr::opts_chunk$set(echo = TRUE)
library(rtpcr)

## ----eval= T------------------------------------------------------------------
# Applying the efficiency function
efficiency(data_efficiency)

## ----eval= T, fig.cap = "relative expression (DDCt) of two different genes. RE tables of any number of genes can be combined and used as input data frame for `plotTwoFactor` function."----

# Example with three factors and one reference gene, without a blocking factor
ANOVA_DCt(data_3factor, numberOfrefGenes = 1, block = NULL)
# Example with a blocking factor
ANOVA_DCt(data_2factorBlock, block = "Block", numberOfrefGenes = 1)

## ----eval= T------------------------------------------------------------------
# Example using two factors with a blocking factor and ANCOVA
ANOVA_DDCt(data_2factorBlock,
numberOfrefGenes = 1,
mainFactor.column = 1,
block = "block",
analysisType = "ancova")

## ----eval= T------------------------------------------------------------------
# Example for repeated measures data (time is the factor of interest)
REPEATED_DDCt(data_repeated_measure_1,
numberOfrefGenes = 1,
factor = "time", block = NULL)

## ----eval= T------------------------------------------------------------------
# Example for unpaired t-test based analysis
TTEST_DDCt(data_ttest,
paired = FALSE,
var.equal = TRUE,
numberOfrefGenes = 1)

## ----eval= T------------------------------------------------------------------
# Before plotting, the result needs to be extracted as below:
res <- ANOVA_DCt(data_1factor, numberOfrefGenes = 1, block = NULL)$Result

# Bar plot
plotOneFactor(res, 1, 2, 7, 8, 11,
    show.groupingLetters = TRUE)

## ----eval= T------------------------------------------------------------------
a <- ANOVA_DCt(data_2factorBlock, block = "Block", numberOfrefGenes = 1)
data <- a$Results

plotTwoFactor(
  data = data,
  x_col = 2,
  y_col = 3,
  group_col = 1,
  Lower.se_col = 8,
  Upper.se_col = 9,
  letters_col = 12,
  show.groupingLetters = TRUE)

## ----eval= T, fig.height = 7, fig.width = 12.5, fig.align = 'center'----------
res <- ANOVA_DCt(data_3factor, 
      numberOfrefGenes = 1, 
      block = NULL)
      
data <- res$Results

plotThreeFactor(data, 
                3,   # x-axis factor
                5,   # bar height
                1,   # fill groups
                2,   # facet grid
                11,   # lower SE column
                12,   # upper SE column
                letters_col = 13,
                show.groupingLetters = TRUE)

## ----eval= T------------------------------------------------------------------

# Returning fold change values of Conc levels sliced by Type
# Returning fold change values from a fitted model.
# Firstly, result of `qpcrANOVAFC` or `qpcrREPEATED` is 
# acquired which includes a model object:
# Assume 'res' is the result from ANOVA_DDCt
res <- ANOVA_DDCt(data_3factor, numberOfrefGenes = 1, mainFactor.column = 1, block = NULL)

Means_DDCt(res$lm_ANOVA, specs = "Conc * Type")

res2 <- Means_DDCt(res$lm_ANOVA, specs = "Conc | Type")

# Returning fold change values of Conc levels sliced by Type*SA interaction
Means_DDCt(res$lm_ANOVA, specs = "Conc | (Type*SA)")

## ----eval= T------------------------------------------------------------------
citation("rtpcr")

