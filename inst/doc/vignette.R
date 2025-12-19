## ----setup, include = FALSE, fig.align='center', warning = F, message=F-------
options(tinytex.verbose = TRUE)
knitr::opts_chunk$set(echo = TRUE)
library(rtpcr)

## ----eval= F------------------------------------------------------------------
# # Installing from CRAN
# install.packages("rtpcr")
# 
# # Loading the package
# library(rtpcr)

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
analysisType = "anova")

## ----eval= T------------------------------------------------------------------
# Example for repeated measures data (time is the factor of interest)
REPEATED_DDCt(data_repeated_measure_1,
numberOfrefGenes = 1,
factor = "time", 
calibratorLevel = "1",
block = NULL)

## ----eval= T------------------------------------------------------------------
# Example for unpaired t-test based analysis

p1 <- TTEST_DDCt(data_1factor_one_ref,
               numberOfrefGenes = 1,
               paired = FALSE,
               var.equal = TRUE,
               plotType = "RE")
p2 <- TTEST_DDCt(data_1factor_one_ref,
                 numberOfrefGenes = 1,
                 paired = FALSE,
                 var.equal = TRUE,
                 plotType = "log2FC")

p1 <- p1$plot
p2 <- p2$plot
rtpcr::multiplot(p1,p2,cols = 2)

## ----eval= T------------------------------------------------------------------
# Before plotting, the result table needs to be extracted:
p1 <- TTEST_DDCt(data_1factor_one_ref,
                 numberOfrefGenes = 1,
                 paired = FALSE,
                 var.equal = TRUE,
                 plotType = "RE")

p2 <- TTEST_DDCt(data_1factor_one_ref,
                 numberOfrefGenes = 1,
                 paired = FALSE,
                 var.equal = TRUE,
                 plotType = "log2FC")

d1 <- p1$Result

# preserve order in plot as in data
d1$Gene <- factor(d1$Gene, levels = unique(d1$Gene))

pl1 <- plotOneFactor(
  d1,
  x_col = 1,
  y_col = 2,
  Lower.se_col = 8,
  Upper.se_col = 9,
  letters_col = 10,
  letters_d = 0.2,
  col_width = 0.8,
  err_width = 0.15,
  fill_colors = "cyan",
  colour = "black",
  alpha = 1,
  base_size = 12,
  legend_position = "none")

pl2 <- plotOneFactor(
  d1,
  x_col = 1,
  y_col = 7,
  Lower.se_col = 11,
  Upper.se_col = 12,
  letters_col = 10,
  letters_d = 0.2,
  col_width = 0.8,
  err_width = 0.15,
  fill_colors = "cyan",
  colour = "black",
  alpha = 1,
  base_size = 12,
  legend_position = "none")

multiplot(pl1, pl2, cols =  2)

## ----eval= T, warning = F, fig.height = 7, fig.width = 12.5, fig.align = 'center', warning = F----
a <- ANOVA_DCt(data_2factorBlock, block = "Block", numberOfrefGenes = 1)
data <- a$Results

p1 <- plotTwoFactor(
  data = data,
  x_col = 2,
  y_col = 3,
  group_col = 1,
  Lower.se_col = 8,
  Upper.se_col = 9,
  letters_col = 12,
  letters_d = 0.2,
  fill_colors = c("aquamarine4", "gold2"),
  alpha = 1,
  col_width = 0.7,
  dodge_width = 0.7,
  base_size = 14, 
  legend_position = c(0.2, 0.8),
  color = "black"
)

library(ggplot2)
p1 <- p1 + scale_y_continuous(breaks = seq(0, 3.6, by = 1),
                        limits = c(0, 3.6),
                        expand = c(0, 0)) + 
  theme(axis.text.x = element_text(size = 14, color = "black", angle = 45),
        axis.text.y = element_text(size = 14,color = "black", angle = 0, hjust = 0.5)) +
  theme(legend.text = element_text(colour = "black", size = 14),
        legend.background = element_rect(fill = "transparent")) +
  ylab("Relative Expression (DCt method)")


p2 <- plotTwoFactor(
  data = data,
  x_col = 2,
  y_col = 4,
  group_col = 1,
  Lower.se_col = 10,
  Upper.se_col = 11,
  letters_col = 12,
  letters_d = 0.2,
  fill_colors = c("aquamarine4", "gold2"),
  alpha = 1,
  col_width = 0.7,
  dodge_width = 0.7,
  base_size = 14, 
  legend_position = c(0.2, 0.8),
  color = "black"
)


p2 <- p2 + 
  scale_y_continuous(expand = c(-1.44, +1.5)) + 
  theme(axis.text.x = element_text(size = 14, color = "black", angle = 45),
        axis.text.y = element_text(size = 14,color = "black", angle = 0, hjust = 0.5)) +
  theme(legend.text = element_text(colour = "black", size = 14),
        legend.background = element_rect(fill = "transparent")) +
  #geom_hline(yintercept = 0, color = "black",  linewidth = 0.6, linetype = "solid") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  ylab("Log2 fold change (DCt method)")

multiplot(p1, p2, cols =  2)

## ----eval= T, fig.height = 7, fig.width = 12.5, fig.align = 'center', warning = F----
res <- ANOVA_DCt(data_3factor, 
      numberOfrefGenes = 1, 
      block = NULL)
      
data <- res$Results

p <- plotThreeFactor(
  data,
  x_col = 3,        # x-axis factor
  y_col = 5,        # bar height
  group_col = 1,    # grouping (fill)
  facet_col = 2,    # faceting factor
  Lower.se_col = 11,
  Upper.se_col = 12,
  letters_col = 13,
  letters_d = 0.3,
  col_width = 0.7, 
  dodge_width = 0.7,# controls spacing
  fill_colors = c("blue", "brown"),
  base_size = 16, 
  alpha = 1,
  legend_position = c(0.1, 0.2)
)

library(ggplot2)
p + theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
)

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

## ----eval= F, warning = F-----------------------------------------------------
# p <- plotOneFactor(...)
# p +
#   geom_hline(yintercept = 0, linetype = "dashed")

## ----eval= F, warning = F-----------------------------------------------------
# p <- plotOneFactor(...)
# p +
#   scale_y_continuous(limits = c(0, 20))

## ----eval= F, warning = F-----------------------------------------------------
# p <- plotTwoFactor(...)
# p +
#   scale_x_discrete(labels = c("A" = "Control", "B" = "Treatment"))

## ----eval= F, warning = F-----------------------------------------------------
# p <- plotTwoFactor(...)
# p +
#   scale_fill_brewer(palette = "Set2")

## ----eval= F, warning = F-----------------------------------------------------
# p <- plotOneFactor(...)
# p +
#   geom_hline(yintercept = 0, linetype = "dashed")

## ----eval= F, warning = F-----------------------------------------------------
# plotOneFactor(...) +
#   geom_hline(yintercept = 0, linetype = "dashed")

## ----eval= F, warning = F-----------------------------------------------------
# Example 1
# p2 <- plotTwoFactor(
#   data = data,
#   x_col = 2,
#   y_col = 4,
#   group_col = 1,
#   Lower.se_col = 10,
#   Upper.se_col = 11,
#   letters_col = 12,
#   letters_d = 0.2,
#   fill_colors = c("aquamarine4", "gold2"),
#   alpha = 1,
#   col_width = 0.7,
#   dodge_width = 0.7,
#   base_size = 16,
#   legend_position = c(0.2, 0.8)
# )
# 
# library(ggplot2)
# p2 + scale_y_continuous(expand = c(-1.5, +1.5)) +
#   theme(axis.text.x = element_text(size = 14, color = "black", angle = 45),
#         axis.text.y = element_text(size = 14,color = "black", angle = 0, hjust = 0.5)) +
#   theme(legend.text = element_text(colour = "black", size = 14),
#         legend.background = element_rect(fill = "transparent"))

## ----eval= T, eval= T, fig.height = 5, fig.width = 10, fig.align = 'center', fig.cap = "QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions."----

residuals <- ANOVA_DCt(data_1factor, numberOfrefGenes = 1, block = NULL)$lmCRD$residuals
shapiro.test(residuals) 

par(mfrow = c(1,2))
plot(residuals)
qqnorm(residuals)
qqline(residuals, col = "red")

## ----eval= T, eval= T, fig.height = 4, fig.width = 4, fig.align = 'center', fig.cap = "QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions."----

a <- REPEATED_DDCt(data_repeated_measure_2, 
                  numberOfrefGenes = 1,
                  calibratorLevel = "1",
                  factor = "time", 
                  block = NULL)


residuals(a$lm)
plot(residuals(a$lm))
qqnorm(residuals(a$lm))
qqline(residuals(a$lm), col = "red")

## ----eval= T------------------------------------------------------------------
# See example input data frame:
data_withTechRep

# Calculating mean of technical replicates
meanTech(data_withTechRep, groups = 1:4)

## ----eval= T, fig.height = 4, fig.width = 5, fig.align = 'center', fig.cap = "Fold change expression of two different genes. FC tables of any number of genes can be combined and used as input data frame for `twoFACTORplot` function."----

a <- REPEATED_DDCt(data_repeated_measure_1,
                   numberOfrefGenes = 1,
                   calibratorLevel = "1",
                   factor = "time", 
                   block = NULL,
                   plot = F)

b <- REPEATED_DDCt(data_repeated_measure_2,
                   factor = "time",
                   calibratorLevel = "1",
                   numberOfrefGenes = 1, 
                   block = NULL,
                   plot = F)


a1 <- a$Relative_Expression_table
b1 <- b$Relative_Expression_table
  
c <- rbind(a1, b1)
c$gene <- factor(c(1,1,1,2,2,2))
c
c$gene <- factor(c$gene, levels = unique(c$gene))
plotTwoFactor(
  c,
  x_col = 1,
  y_col = 2,
  group_col = 13,
  Lower.se_col = 9,
  Upper.se_col = 10,
  letters_col = 5,
  letters_d = 0.2,
  dodge_width = 0.7,
  col_width = 0.7,
  err_width = 0.15,
  fill_colors = c("aquamarine4", "gold2"),,
  alpha = 1,
  legend_position = c(0.2, 0.8),
  base_size = 12
)

## ----eval= T------------------------------------------------------------------
citation("rtpcr")

