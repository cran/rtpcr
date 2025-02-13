## ----setup, include = FALSE, fig.align='center', warning = F, message=F-------
options(tinytex.verbose = TRUE)
knitr::opts_chunk$set(echo = TRUE)

## ----eval= T, include= F, message=FALSE, warning = FALSE----------------------
library(rtpcr)
library(multcomp)
library(dplyr)
library(reshape2)
library(tidyr)
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

## ----eval = T, , fig.height = 3, fig.width = 5, fig.align = 'center', fig.cap = 'Standard curve and the amplification efficiency analysis of genes. Required iput data include dilutions and Ct value columns for different genes.', warning = FALSE, message = FALSE----
efficiency(data_efficiency)

## ----eval= T, fig.height = 3, fig.width = 5, fig.align = 'center'-------------
data_ttest

## ----eval= T------------------------------------------------------------------
qpcrTTEST(data_ttest, 
          numberOfrefGenes = 1,
          paired = F, 
          var.equal = T)

## ----eval= T, fig.height=3, fig.width=8, fig.align='center', fig.cap = "Average Fold changes of three target genes relative to the control condition computed by unpaired t-tests via `qpcrTTESTplot` function. Confidence interval (ci) and standard error (se) has been used as error bar in 'A' and 'B', respectively.", warning = F, message = F----

# Producing the plot
t1 <- qpcrTTESTplot(data_ttest,
              numberOfrefGenes = 1,
              fontsizePvalue = 4,
              errorbar = "ci")

# Producing the plot: specifying gene order
t2 <- qpcrTTESTplot(data_ttest,
              numberOfrefGenes = 1,
              order = c("C2H2-01", "C2H2-12", "C2H2-26"),
              paired = FALSE,
              var.equal = TRUE,
              width = 0.5,
              fill = "palegreen",
              y.axis.adjust = 3,
              y.axis.by = 2,
              ylab = "Average Fold Change",
              xlab = "none",
              fontsizePvalue = 4)

multiplot(t1, t2, cols = 2)
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))

## ----eval = T, fig.height = 3, fig.width = 5, fig.align='center', fig.cap = "Statistical table and figure of the Fold change expression of a gene in three different levels of Drough stress relative to the D0 as reference or calibrator level produced by the `qpcrANOVAFC` function. The other factor i.e. Genotype has been concidered as covariate."----
# See sample data
data_2factor

qpcrANOVAFC(data_2factor, 
           numberOfrefGenes = 1, 
           block = NULL,
           analysisType = "ancova",
           mainFactor.column = 2,
           fontsizePvalue = 4,
           x.axis.labels.rename = "none")

## ----eval= T------------------------------------------------------------------
# See a sample dataset
data_3factor

## ----eval= T, fig.height = 3, fig.width = 5-----------------------------------
# If the data include technical replicates, means of technical replicates
# should be calculated first using meanTech function.

# Applying ANOVA analysis
res <- qpcrANOVARE(data_2factor,
                   numberOfrefGenes = 1,
                   block = NULL)
res$Result
res$Post_hoc_Test

## ----eval= T, fig.height = 4, fig.width = 9, fig.align = 'center', fig.cap = "A: bar plot representing Relative expression of a gene under three levels of a factor generated using `oneFACTORplot` function, B: Plot of the Fold change expression produced by the `qpcrANOVAFC` function from the same data used for 'A'. The first element in the `mainFactor.level.order` argument (here L1) is served as the Reference level, although the x-axis names have later been renamed by the `x.axis.labels.rename` argument. Error bars represent 95% confidence interval in A and standard error in B."----

# Before plotting, the statistical analysis should be done:
out2 <- qpcrANOVARE(data_1factor, numberOfrefGenes = 1, block = NULL)$Result

f1 <- oneFACTORplot(out2,
              width = 0.2,
              fill = "skyblue",
              y.axis.adjust = 0.5,
              y.axis.by = 1,
              errorbar = "ci",
              show.letters = TRUE,
              letter.position.adjust = 0.1,
              ylab = "Relative Expression",
              xlab = "Factor Levels",
              fontsize = 12,
              fontsizePvalue = 4)


addline_format <- function(x,...){
    gsub('\\s','\n',x)
}
f2 <- qpcrANOVAFC(data_1factor,
                 numberOfrefGenes = 1,
                 mainFactor.column = 1,
                 block = NULL,
                 mainFactor.level.order = c("L1","L2","L3"),
                 width = 0.5,
                 fill = c("skyblue","#79CDCD"),
                 y.axis.by = 1,
                 letter.position.adjust = 0,
                 y.axis.adjust = 1,
                 ylab = "Fold Change",
                 fontsize = 12, plot = F,
                 x.axis.labels.rename = addline_format(c("Control", 
                                                       "Treatment_1 vs Control", 
                                                       "Treatment_2 vs Control")))


multiplot(f1, f2$FC_Plot_of_the_main_factor_levels, cols = 2)
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))

## ----eval= T, include = T, fig.height = 4, fig.width = 9, fig.align = 'center', fig.cap = "Relative expression of a target gene under two different factors of genotype (with two levels) and drought (with three levels). Error bars represent standard error. Means (columns) lacking letters in common have significant difference at alpha = 0.05 as resulted from a `LSD.test`."----

# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVARE(data_2factor, numberOfrefGenes = 1, block = NULL)$Result
Final_data <- qpcrANOVARE(data_2factor, numberOfrefGenes = 1, block = NULL)$Final_data

# Plot of the 'res' data with 'Genotype' as grouping factor
q1 <- twoFACTORplot(res,
   x.axis.factor = Drought,
   group.factor = Genotype,
   errorbar = "se",
   width = 0.5,
   fill = "Greens",
   y.axis.adjust = 0.5,
   y.axis.by = 2,
   ylab = "Relative Expression",
   xlab = "Drought Levels",
   legend.position = c(0.15, 0.8),
   show.letters = TRUE,
   fontsizePvalue = 4)

# Plotting the same data with 'Drought' as grouping factor
q2 <- twoFACTORplot(res,
   x.axis.factor = Genotype,
   group.factor = Drought,
   errorbar = "se",
   xlab = "Genotype",
   fill = "Blues",
   legend.position = c(0.15, 0.8),
   show.letters = FALSE,
   show.errorbars = T,
   fontsizePvalue = 4)

multiplot(q1, q2, cols = 2)
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))

## ----fig.height = 5, fig.width = 11, fig.align = 'center', fig.cap = "A and B) Relative expression (RE) of a target gene from a three-factorial experiment data produced by  `threeFACTORplot`  function. Error bars represent standard error (A), although can be set to confidence interval (B). Means (columns) lacking letters in common have significant differences at alpha = 0.05 as resulted from an ‘LSD.test’."----
# Before plotting, the result needs to be extracted as below:
res <- qpcrANOVARE(data_3factor, numberOfrefGenes = 1, block = NULL)$Result
res

# releveling a factor levels first
res$Conc <- factor(res$Conc, levels = c("L","M","H"))
res$Type <- factor(res$Type, levels = c("S","R"))

# Arrange the first three colunms of the result table.
# This determines the columns order and shapes the plot output.
p1 <- threeFACTORplot(res,
    arrangement = c(3, 1, 2),
    errorbar = "se",
    legend.position = c(0.2, 0.85),
    xlab = "condition",
    fontsizePvalue = 4)


# When using ci as error, increase y.axis.adjust to see the plot correctly!
p2 <- threeFACTORplot(res,
   arrangement = c(2, 3, 1),
   bar.width = 0.8,
   fill = "Greens",
   xlab = "Drought",
   ylab = "Relative Expression",
   errorbar = "ci",
   y.axis.adjust = 2,
   y.axis.by = 2,
   letter.position.adjust = 0.6,
   legend.title = "Genotype",
   fontsize = 12,
   legend.position = c(0.2, 0.8),
   show.letters = TRUE,
   fontsizePvalue = 4)

multiplot(p1, p2, cols = 2)
grid.text("A", x = 0.02, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))
grid.text("B", x = 0.52, y = 1, just = c("right", "top"), gp=gpar(fontsize=16))

## ----eval=T, fig.height = 4, fig.width = 7, fig.align = 'center', fig.cap = "Fold change expression (FC) of a target gene from a one and a two factorial experiment data produced by  `qpcrREPEATED`  function. Error bars represent standard error (A), although can be set to confidence interval."----

a <- qpcrREPEATED(data_repeated_measure_1,
                  numberOfrefGenes = 1,
                  block = NULL,
                  fill = c("#778899", "#BCD2EE"),
                  factor = "time",
                  axis.text.x.angle = 45,
                  axis.text.x.hjust = 1,
                  plot = F)

b <- qpcrREPEATED(data_repeated_measure_2,
                  numberOfrefGenes = 1,
                  factor = "time",
                  block = NULL,
                  axis.text.x.angle = 45,
                  axis.text.x.hjust = 1,
                  plot = F)

multiplot(a, b, cols = 2)

## ----eval=T-------------------------------------------------------------------
# Returning fold change values from a fitted model.
# Firstly, result of `qpcrANOVAFC` or `qpcrREPEATED` is 
# acquired which includes a model object:
res <- qpcrANOVAFC(data_3factor, numberOfrefGenes = 1, mainFactor.column = 1, block = NULL)

# Returning fold change values of Conc levels from a fitted model:
qpcrMeans(res$lm_ANOVA, specs = "Conc")


# Returning fold change values of Conc levels sliced by Type*SA:
qpcrMeans(res$lm_ANOVA, specs = "Conc | (Type*SA)")

# Returning fold change values of Conc
qpcrMeans(res$lm_ANOVA, specs = "Conc * Type")

# Returning fold change values of Conc levels sliced by Type:
res2 <- qpcrMeans(res$lm_ANOVA, specs = "Conc | Type")
twoFACTORplot(res2, x.axis.factor = contrast, ylab = "Fold Change",
              group.factor = Type, errorbar = "ci")

## ----eval=T, include = T, fig.height = 4, fig.width = 6, fig.align = 'center'----

library(ggplot2)
b <- qpcrANOVARE(data_3factor, numberOfrefGenes = 1, block = NULL)$Result
a <- qpcrANOVARE(data_3factor, numberOfrefGenes = 1, block = NULL)$Final_data

# Arrange factor levels to your desired order:
b$Conc <- factor(b$Conc, levels = c("L","M","H"))
a$Conc <- factor(a$Conc, levels = c("L","M","H"))

# Generating plot
ggplot(b, aes(x = Type, y = RE, fill = factor(Conc))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ SA) +
  scale_fill_brewer(palette = "Reds") +
  xlab("Type") +
  ylab("Relative Expression") +
  geom_point(data = a, aes(x = Type, y = (2^(-wDCt)), fill = factor(Conc)), 
             position = position_dodge(width = 0.9), color = "black") +
  ylab("ylab") +
  xlab("xlab") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
        axis.title  = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  theme(legend.position  = c(0.2, 0.7)) +
  theme(legend.title = element_text(size = 12, color = "black")) +
  scale_y_continuous(breaks = seq(0, max(b$RE) + max(b$se) + 0.1, by = 5), 
                     limits = c(0, max(b$RE) + max(b$se) + 0.1), expand = c(0, 0)) 

## ----eval= T, eval= T, fig.height = 5, fig.width = 10, fig.align = 'center', fig.cap = "QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions."----

residuals <- qpcrANOVARE(data_1factor, numberOfrefGenes = 1, block = NULL)$lmCRD$residuals
shapiro.test(residuals) 

par(mfrow = c(1,2))
plot(residuals)
qqnorm(residuals)
qqline(residuals, col = "red")

## ----eval= T, eval= T, fig.height = 4, fig.width = 4, fig.align = 'center', fig.cap = "QQ-plot for the normality assessment of the residuals derived from `t.test` or `lm` functions."----

a <- qpcrREPEATED(data_repeated_measure_2, 
                  numberOfrefGenes = 1, 
                  factor = "time", 
                  block = NULL,
                  y.axis.adjust = 1.5)


residuals(a$lm)
plot(residuals(a$lm))
qqnorm(residuals(a$lm))
qqline(residuals(a$lm), col = "red")

## ----eval= T------------------------------------------------------------------
# See example input data frame:
data_withTechRep

# Calculating mean of technical replicates
meanTech(data_withTechRep, groups = 1:4)

## ----eval= T, eval= T, , fig.height = 4, fig.width = 5, fig.align = 'center', fig.cap = "Fold change expression of two different genes. FC tables of any number of genes can be combined and used as input data frame for `twoFACTORplot` function."----

a <- qpcrREPEATED(data_repeated_measure_1,
             numberOfrefGenes = 1,
             factor = "time", block = NULL)

b <- qpcrREPEATED(data_repeated_measure_2,
                  factor = "time",
                  numberOfrefGenes = 1, block = NULL)


a1 <- a$FC_statistics_of_the_main_factor
b1 <- b$FC_statistics_of_the_main_factor

c <- rbind(a1, b1)
c$gene <- factor(c(1,1,1,2,2,2))
c

twoFACTORplot(c, x.axis.factor = contrast, 
              group.factor = gene, fill = 'Reds', errorbar = "se",
              ylab = "FC", axis.text.x.angle = 45, y.axis.adjust = 1.5,
              axis.text.x.hjust = 1, legend.position = c(0.2, 0.8))

## ----eval= F, include = T, fig.height = 4, fig.width = 5----------------------
#  
#  b <- qpcrANOVAFC(data_2factor,
#              numberOfrefGenes = 1,
#              mainFactor.column = 1,
#              block = NULL,
#              mainFactor.level.order = c("S", "R"),
#              fill = c("#CDC673", "#EEDD82"),
#              analysisType = "ancova",
#              fontsizePvalue = 7,
#              y.axis.adjust = 0.1, width = 0.35)
#  
#  
#  
#  library(ggplot2)
#  p2 <- b$FC_Plot_of_the_main_factor_levels
#  p2 + theme_set(theme_classic(base_size = 20))
#  

## ----eval= F------------------------------------------------------------------
#  citation("rtpcr")

