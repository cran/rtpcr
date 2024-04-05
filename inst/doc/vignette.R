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
library(grid)

## ----eval= F, include= T, message=FALSE, warning = FALSE----------------------
#  
#  # install `rtpcr` from github (under development)
#  
#  devtools::install_github("mirzaghaderi/rtpcr")
#  
#  # I strongly recommend to install the package with the vignette as it contains information about how to use the 'rtpcr' package. Through the following code, Vignette is installed as well.
#  
#  devtools::install_github("mirzaghaderi/rtpcr", build_vignettes = TRUE)

## ----eval= T------------------------------------------------------------------
data_efficiency

## ----eval = T, fig.height = 3, fig.align = 'center', fig.cap = 'Standard curve and the amplification efficiency analysis of target and reference genes. A sample data arrangement that is required as input for the calculation of amplification efficiency by the efficiency function.', warning = FALSE, message = FALSE----
efficiency(data_efficiency)

## ----eval= T, fig.height = 3, fig.width = 5, fig.align = 'center'-------------
data_ttest

## ----eval= T------------------------------------------------------------------
qpcrTTEST(data_ttest, 
          numberOfrefGenes = 1,
          paired = F, 
          var.equal = T)

## ----eval= T, fig.height=3, fig.width=8, fig.align='center', fig.cap = "Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests via ‘qpcrTTESTplot’ function.", warning = F, message = F----

# Producing the plot
t1 <- qpcrTTESTplot(data_ttest,
              numberOfrefGenes = 1)

# Producing the plot: specifying gene order
t2 <- qpcrTTESTplot(data_ttest,
              numberOfrefGenes = 1,
              order = c("C2H2-01", "C2H2-12", "C2H2-26"),
              paired = FALSE,
              var.equal = TRUE,
              width = 0.5,
              fill = "palegreen",
              y.axis.adjust = 0,
              y.axis.by = 2,
              ylab = "Average Fold Change (FC)",
              xlab = "Gene")

multiplot(t1, t2, cols = 2)
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))

## ----eval = T, fig.height = 3, fig.width = 5, fig.align='center'--------------
# See sample data
data_2factor

qpcrANCOVA(data_2factor, 
           numberOfrefGenes = 1, 
           analysisType = "ancova", 
           main.factor = 2, 
           y.axis.adjust = 0.3,
           levels = c(3, 2, 1))

## ----eval= T------------------------------------------------------------------
# See a sample dataset
data_3factor_a

## ----eval= T, fig.height = 3, fig.width = 5-----------------------------------
# If the data include technical replicates, means of technical replicates
# should be calculated first using meanTech function.

# Applying ANOVA analysis
res <- qpcrANOVA(data_2factor,
          numberOfrefGenes = 1,
          p.adj = "none")
res$Result
res$Post_hoc_Test

## ----eval= T, fig.height = 4, fig.width = 9, fig.align = 'center', fig.cap = "A) A bar plot representing Relative expression of a gene under three levels of a factor generated using ‘oneFACTORplot’ function, B) Plot of average Fold changes produced by the ‘qpcrANCOVA’ function from the same data as ‘C’. Check level can be changed by user. Error bars represent 95% confidence interval."----

# Before plotting, the result needs to be extracted as below:
out2 <- qpcrANOVA(data_1factor, numberOfrefGenes = 1)$Result

f1 <- oneFACTORplot(out2,
              width = 0.2,
              fill = "skyblue",
              y.axis.adjust = 0.5,
              y.axis.by = 1,
              errorbar = "ci",
              show.letters = TRUE,
              letter.position.adjust = 0.1,
              ylab = "Relative Expression (RE)",
              xlab = "Factor Levels",
              fontsize = 12)

f2 <- qpcrANCOVA(data_1factor, 
           numberOfrefGenes = 1,
           analysisType = "ancova", 
           main.factor = 1,
           levels = c(1, 2, 3))

multiplot(f1, f2, cols = 2)
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))

## ----eval= T, include = T, fig.height = 4, fig.width = 9, fig.align = 'center', fig.cap = "Average relative expression of a target gene under two different factors of genotype (with two levels) and drought (with three levels). Error bars represent standard deviations. Means (columns) lacking letters in common have significant difference at alpha = 0.05 as resulted from the `LSD.test` of agricolae package."----

# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVA(data_2factor, numberOfrefGenes = 1)

# Plot of the 'res' data with 'Genotype' as grouping factor
q1 <- twoFACTORplot(res,
   x.axis.factor = Drought,
   group.factor = Genotype,
   width = 0.5,
   fill = "Greens",
   y.axis.adjust = 0.5,
   y.axis.by = 2,
   ylab = "Relative Expression",
   xlab = "Drought Levels",
   legend.position = c(0.15, 0.8),
   show.letters = TRUE)

# Plotting the same data with 'Drought' as grouping factor
q2 <- twoFACTORplot(res,
   x.axis.factor = Genotype,
   group.factor = Drought,
   xlab = "Genotype",
   fill = "Blues",
   legend.position = c(0.15, 0.8),
   show.letters = FALSE,
   show.errorbars = F,
   show.points = T)

multiplot(q1, q2, cols = 2)
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))

## ----fig.height = 5, fig.width = 11, fig.align = 'center', fig.cap = "A and B) Relative expression (RE) of a target gene under two or three factors produced by ‘twoFACTORplot’ and ‘threeFACTORplot’ functions, respectively. Error bars represent standard deviations (can be set to confidence interval). Means (columns) lacking letters in common have significant differences at alpha = 0.05 as resulted from an ‘LSD.test’."----
# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVA(data_3factor_b, numberOfrefGenes = 1)$Result
res

# releveling a factor levels first
res$Conc <- factor(res$Conc, levels = c("L","M","H"))
res$Type <- factor(res$Type, levels = c("S","R"))

# Arrange the first three colunms of the result table.
# This determines the columns order and shapes the plot output.
p1 <- threeFACTORplot(res,
    arrangement = c(3, 1, 2),
    legend.position = c(0.2, 0.85),
    xlab = "condition")


# When using ci as error, increase y.axis.adjust to see the plot correctly!
p2 <- threeFACTORplot(res,
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

multiplot(p1, p2, cols = 2)
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))

## ----eval=F, include = T------------------------------------------------------
#  
#  b <- qpcrANOVA(data_3factor_a, numberOfrefGenes = 1)$Result
#  a <- qpcrANOVA(data_3factor_a, numberOfrefGenes = 1)$Final_data
#  
#  ggplot(b, aes(x = Genotype, y = RE, fill = factor(Drought))) +
#    geom_bar(stat = "identity", position = "dodge") +
#    facet_wrap(~ SA) +
#    scale_fill_brewer(palette = "Reds") +
#    xlab("Genotype") +
#    ylab("Relative Expression") +
#    geom_point(data = a, aes(x = Genotype, y = (10^(-wDCt)), fill = factor(Drought)),
#               position = position_dodge(width = 0.9), color = "black") +
#    ylab("ylab") +
#    xlab("xlab") +
#    theme_bw() +
#    theme(axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
#          axis.text.y = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
#          axis.title  = element_text(size = 12),
#          legend.text = element_text(size = 12)) +
#    theme(legend.position  = c(0.2, 0.7)) +
#    theme(legend.title = element_text(size = 12, color = "black")) +
#    scale_y_continuous(breaks = seq(0, max(b$RE) + max(b$std) + 0.1, by = 5),
#                       limits = c(0, max(b$RE) + max(b$std) + 0.1), expand = c(0, 0))

## ----eval= T, eval= T, , fig.height = 4, fig.width = 5, fig.align = 'center', fig.cap = "QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions."----

residuals <- qpcrANOVA(data_1factor, numberOfrefGenes = 1)$lmCRD$residuals
shapiro.test(residuals) 

qqnorm(residuals)
qqline(residuals, col = "red")

## ----eval= T, eval= T---------------------------------------------------------
# See example input data frame:
data_withTechRep

# Calculating mean of technical replicates
meanTech(data_withTechRep, groups = 1:4)

## ----eval= F------------------------------------------------------------------
#  citation("rtpcr")

