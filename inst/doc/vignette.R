## ----setup, include = FALSE, fig.align='center', warning = F, message=F-------
options(tinytex.verbose = TRUE)
knitr::opts_chunk$set(echo = TRUE)

## ----eval= T, include= F, message=FALSE, warning = FALSE----------------------
library(rtpcr)
library(agricolae)
library(dplyr)
library(reshape2)
library(tidyr)
library(lme4)
library(ggplot2)

## ----eval= F, include= T, message=FALSE, warning = FALSE----------------------
#  install.packages("rtpcr")
#  library(rtpcr)

## ----eval= F, include= T, message=FALSE, warning = FALSE----------------------
#  # install `rtpcr` from github (under development)
#  devtools::install_github("mirzaghaderi/rtpcr", build_vignettes = TRUE)

## ----eval= T------------------------------------------------------------------
data_efficiency

## ----eval = T, fig.height = 3, fig.align = 'center', fig.cap = 'Standard curve and the amplification efficiency analysis of target and reference genes.', warning = FALSE, message = FALSE----
efficiency(data_efficiency)

## ----eval= T, fig.height = 3, fig.width = 5, fig.align = 'center'-------------
data_ttest

## ----eval= T------------------------------------------------------------------
qpcrTTEST(data_ttest, 
          numberOfrefGenes = 1,
          paired = F, 
          var.equal = T)

## ----eval= T, fig.height=3, fig.width=5, fig.align='center', fig.cap = "Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests. Error bars represent 95% confidence interval.", warning = F, message = F----

# Producing the plot
qpcrTTESTplot(data_ttest,
              numberOfrefGenes = 1)

# Producing the plot: specifying gene order
qpcrTTESTplot(data_ttest,
              numberOfrefGenes = 1,
              order = c("C2H2-01", "C2H2-12", "C2H2-26"),
              paired = FALSE,
              var.equal = TRUE,
              width = 0.5,
              fill = "skyblue",
              y.axis.adjust = 0,
              y.axis.by = 2,
              ylab = "Average Fold Change",
              xlab = "Gene")

## ----eval= T------------------------------------------------------------------
# See a sample dataset
data_3factor_a

## ----eval= T, fig.height = 3, fig.width = 5-----------------------------------
# If the data include technical replicates, means of technical replicates
# should be calculated first using meanTech function.

# Applying ANOVA analysis
qpcrANOVA(data_3factor_a,
          numberOfrefGenes = 1,
          p.adj = "none")

## ----eval= T, fig.height = 4, fig.width = 6, fig.align = 'center', fig.cap = "Average relative expression of a target gene under two different factors of genotype (with two levels) and drought (with three levels). Error bars represent standard deviations. Means (columns) lacking letters in common have significant difference at alpha = 0.05 as resulted from the `LSD.test` of agricolae package."----
# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVA(data_2factor, numberOfrefGenes = 1)$Result
res

# Plot of the 'res' data with 'Genotype' as grouping factor
twoFACTORplot(res,
   x.axis.factor = Drought,
   group.factor = Genotype,
   width = 0.5,
   fill = "Greens",
   y.axis.adjust = 0.5,
   y.axis.by = 2,
   ylab = "Relative Expression",
   xlab = "Drought Levels",
   legend.position = c(0.09, 0.8),
   show.letters = TRUE)

# Plotting the same data with 'Drought' as grouping factor
twoFACTORplot(res,
   x.axis.factor = Genotype,
   group.factor = Drought,
   xlab = "Genotype",
   fill = "Blues",
   show.letters = FALSE)

## ----fig.height = 5, fig.width = 7, fig.align = 'center', fig.cap = "Average relative expression of a target gene under three different conditions with two, two and three levels. Error bars can be standard deviation or confidence interval. Means lacking letters in common have significant difference at alpha = 0.05 resulted from the `LSD.test` of agricolae package."----
# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVA(data_3factor_b, numberOfrefGenes = 1)$Result
res

# Arrange the first three colunms of the result table.
# This determines the columns order and shapes the plot output.
threeFACTORplot(res,
    arrangement = c(3, 1, 2),
    legend.position = c(0.9, 0.85),
    xlab = "condition")



threeFACTORplot(res,
   arrangement = c(1, 2, 3),
   bar.width = 0.5,
   fill = "Greys",
   xlab = "Genotype",
   ylab = "Relative Expression")


# releveling a factor levels first
res$Conc <- factor(res$Conc, levels = c("L","M","H"))
res$Type <- factor(res$Type, levels = c("S","R"))

# Producing the plot
threeFACTORplot(res,
   arrangement = c(2, 3, 1),
   bar.width = 0.5,
   fill = "Reds",
   xlab = "Drought",
   ylab = "Relative Expression",
   errorbar = "std",
   legend.title = "Genotype",
   legend.position = c(0.2, 0.8))


# When using ci as error, increase y.axis.adjust to see the plot correctly!
threeFACTORplot(res,
   arrangement = c(2, 3, 1),
   bar.width = 0.8,
   fill = "Greens",
   xlab = "Drought",
   ylab = "Relative Expression",
   errorbar = "ci",
   y.axis.adjust = 8,
   y.axis.by = 2,
   letter.position.adjust = 0.6,
   legend.title = "Genotype",
   fontsize = 12,
   legend.position = c(0.2, 0.8),
   show.letters = TRUE)

## ----eval= T, eval= T, , fig.height = 4, fig.width = 5, fig.align = 'center', fig.cap = "QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions."----

residualsCRD <- qpcrANOVA(data_3factor_b, numberOfrefGenes = 1)$lmCRD$residuals
shapiro.test(residualsCRD) 
qqnorm(residualsCRD)
qqline(residualsCRD, col = "red")

## ----eval= T, eval= T---------------------------------------------------------
# See example input data frame:
data_withTechRep

# Calculating mean of technical replicates
meanTech(data_withTechRep, groups = 1:4)

## ----eval= F------------------------------------------------------------------
#  citation("rtpcr")

