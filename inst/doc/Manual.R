## ----setup, include = FALSE, fig.align='center', warning = F, message=F-------
options(tinytex.verbose = TRUE)
knitr::opts_chunk$set(echo = TRUE)
library(rtpcr)

## ----eval = F-----------------------------------------------------------------
# # Installing from CRAN
# install.packages("rtpcr")
# 
# # Loading the package
# library(rtpcr)

## ----eval = F-----------------------------------------------------------------
# devtools::install_github("mirzaghaderi/rtpcr", build_vignettes = TRUE)

## ----eval= F------------------------------------------------------------------
# # Applying the efficiency function
# data <- read.csv(system.file("extdata", "data_efficiency1.csv", package = "rtpcr"))
# data
# dilutions	Gene1	Gene2	Gene3
# 1.00	25.58	24.25	22.61
# 1.00	25.54	24.13	22.68
# 1.00	25.50	24.04	22.63
# 0.50	26.71	25.56	23.67
# 0.50	26.73	25.43	23.65
# 0.50	26.87	26.01	23.70
# 0.20	28.17	27.37	25.11
# 0.20	28.07	26.94	25.12
# 0.20	28.11	27.14	25.11
# 0.10	29.20	28.05	26.17
# 0.10	29.49	28.89	26.15
# 0.10	29.07	28.32	26.15
# 0.05	30.17	29.50	27.12
# 0.05	30.14	29.93	27.14
# 0.05	30.12	29.71	27.16
# 0.02	31.35	30.69	28.52
# 0.02	31.35	30.54	28.57
# 0.02	31.35	30.04	28.53
# 0.01	32.55	31.12	29.49
# 0.01	32.45	31.29	29.48
# 0.01	32.28	31.15	29.26
# 
# # Analysis
# efficiency(data)
# 
# $Efficiency
#    Gene     Slope        R2        E
# 1 Gene1 -3.388094 0.9965504 1.973110
# 2 Gene2 -3.528125 0.9713914 1.920599
# 3 Gene3 -3.414551 0.9990278 1.962747
# 
# $Slope_compare
# $contrasts
#  contrast          estimate    SE df t.ratio p.value
#  C2H2.26 - C2H2.01   0.1400 0.121 57   1.157  0.4837
#  C2H2.26 - GAPDH     0.0265 0.121 57   0.219  0.9740
#  C2H2.01 - GAPDH    -0.1136 0.121 57  -0.938  0.6186

## ----eval= F------------------------------------------------------------------
# data1 <- read.csv(system.file("extdata", "data_Yuan2006PMCBioinf.csv", package = "rtpcr"))
# data1
#       Con r target target_Ct Actin Actin_Ct
#   control 1   1.88  21.13000  1.67 20.70333
#   control 2   1.88  21.55667  1.67 20.35000
#   control 3   1.88  21.33000  1.67 20.75333
# treatment 1   1.88  22.09000  1.67 20.24333
# treatment 2   1.88  22.69667  1.67 20.54000
# treatment 3   1.88  23.05333  1.67 20.50000
# 
# # Anova analysis
# ANOVA_DDCt(
#   data,
#   specs = "condition",
#   numOfFactors = 1,
#   numberOfrefGenes = 1,
#   block = NULL)
# 
# 
# # An example of a properly arranged dataset from a repeated-measures experiment.
# data2 <- read.csv(system.file("extdata", "data_repeated_measure_1.csv", package = "rtpcr"))
# data2
# time	id  E_Target	Ct_target   E_Ref      Ct_Ref
#    1	 1	   	2	    18.92	2	32.77
#    1	 2	   	2	    15.82	2	32.45
#    1	 3	   	2	    19.84	2	31.62
#    2	 1	  	2	    19.46	2	33.03
#    2	 2	   	2	    17.56	2	33.24
#    2	 3	   	2	    19.74	2	32.08
#    3	 1	  	2	    15.73	2	32.95
#    3	 2	  	2	    17.21	2	33.64
#    3	 3	   	2	    18.09	2	33.40
# 
# # Repeated measure analysis
# res <- ANOVA_DDCt(
#   data,
#   numOfFactors = 1,
#   numberOfrefGenes = 1,
#   specs = "time",
#   block = NULL, model = wDCt ~ time + (1 | id))
# 
# 
# # Paired t.test (equivalent to repeated measure analysis, but not always the same results, due to different calculation methods!)
# TTEST_DDCt(
#   data2[1:6,],
#   numberOfrefGenes = 1,
#   paired = T)
# 
# 
# # Anova analysis
# data3 <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
# 
# res <- ANOVA_DDCt(
#   x = data3,
#   specs = "Type | Concentration",
#   numOfFactors = 2,
#   numberOfrefGenes = 3,
#   block = "block",
#   analyseAllTarget = TRUE)

## ----eval= F------------------------------------------------------------------
# # Relative expression table for the specified column in the input data:
# data3 <- read.csv(system.file("extdata", "data_2factorBlock3ref.csv", package = "rtpcr"))
# 
# res <- ANOVA_DDCt(
#   x = data3,
#   specs = "Concentration",
#   numOfFactors = 2,
#   numberOfrefGenes = 3,
#   block = "block",
#   analyseAllTarget = TRUE)
# 
# # Relative Expression
# #   gene contrast     ddCt      RE   log2FC     LCL     UCL      se Lower.se.RE Upper.se.RE Lower.se.log2FC Upper.se.log2FC  pvalue sig
# # 1   PO       L1  0.00000 1.00000  0.00000 0.00000 0.00000 0.13940     0.90790     1.10144         0.00000         0.00000 1.00000
# # 2   PO L2 vs L1 -0.94610 1.92666  0.94610 1.25860 2.94934 0.14499     1.74245     2.13036         0.85564         1.04613 0.00116  **
# # 3   PO L3 vs L1 -2.19198 4.56931  2.19198 3.08069 6.77724 0.29402     3.72685     5.60221         1.78783         2.68748 0.00000 ***
# # 4  NLM       L1  0.00000 1.00000  0.00000 0.00000 0.00000 0.91809     0.52921     1.88962         0.00000         0.00000 1.00000
# # 5  NLM L2 vs L1  0.86568 0.54879 -0.86568 0.39830 0.75614 0.36616     0.42577     0.70734        -1.11579        -0.67163 0.00018 ***
# # 6  NLM L3 vs L1 -1.44341 2.71964  1.44341 1.94670 3.79946 0.17132     2.41511     3.06256         1.28179         1.62542 0.00000 ***
# #
# # The L1 level was used as calibrator.
# # Note: Using default model for statistical analysis: wDCt ~ block + Concentration * Type
# 
# 
# ANOVA_table <- res$perGene$PO$ANOVA_table
# ANOVA_table
# 
# lm <- res$perGene$PO$lm
# lm
# 
# lm_formula <- res$perGene$gene_name$lm_formula
# lm_formula
# 
# residuals <- resid(res$perGene$gene_name$lm)
# residuals

## ----eval= F, warning = F, fig.height = 7, fig.width = 12.5, fig.align = 'center', warning = F----
# data <- read.csv(system.file("extdata", "data_3factor.csv", package = "rtpcr"))
# #Perform analysis first
# res <- ANOVA_DCt(
#   data,
#   numOfFactors = 3,
#   numberOfrefGenes = 1,
#   block = NULL)
# 
# df <- res$relativeExpression
# df
# # Generate three-factor bar plot
# plotFactor(
#   df,
#   x_col = "SA",
#   y_col = "log2FC",
#   group_col = "Type",
#   facet_col = "Conc",
#   Lower.se_col = "Lower.se.log2FC",
#   Upper.se_col = "Upper.se.log2FC",
#   letters_col = "sig",
#   letters_d = 0.3,
#   col_width = 0.7,
#   dodge_width = 0.7,
#   fill_colors = c("palegreen3", "skyblue"),
#   color = "black",
#   base_size = 14,
#   alpha = 1,
#   legend_position = c(0.1, 0.2))

## ----eval= F, fig.height = 7, fig.width = 12.5, fig.align = 'center', warning = F----
# data <- read.csv(system.file("extdata", "data_2factorBlock.csv", package = "rtpcr"))
# res <- ANOVA_DCt(data,
#       numOfFactors = 2,
#       block = "block",
#       numberOfrefGenes = 1)
# 
# df <- res$relativeExpression
# 
# plotFactor(
#   data = df,
#   x_col = "factor2",
#   y_col = "RE",
#   group_col = "factor1",
#   Lower.se_col = "Lower.se.RE",
#   Upper.se_col = "Upper.se.RE",
#   letters_col = "sig",
#   letters_d = 0.2,
#   fill_colors = c("aquamarine4", "gold2"),
#   color = "black",
#   alpha = 1,
#   col_width = 0.7,
#   dodge_width = 0.7,
#   base_size = 16,
#   legend_position = c(0.8, 0.8))

## ----eval= F, warning = F-----------------------------------------------------
# # Heffer et al., 2020, PlosOne
# library(dplyr)
# df <- read.csv(system.file("extdata", "data_Heffer2020PlosOne.csv", package = "rtpcr"))
# 
# res <- ANOVA_DDCt(
#   df,
#   numOfFactors = 1,
#   specs = "Treatment",
#   numberOfrefGenes = 1,
#   block = NULL)
# 
# data <- res$relativeExpression
# 
# # Selecting only the first words in 'contrast' column to be used as the x-axis labels.
# data$contrast <- sub(" .*", "", data$contrast)
# 
# plotFactor(
#   data = data,
#   x_col = "contrast",
#   y_col = "RE",
#   group_col = "contrast",
#   facet_col = "gene",
#   Lower.se_col = "Lower.se.RE",
#   Upper.se_col = "Upper.se.RE",
#   letters_col = "sig",
#   legend_position = "none") # more controlling arguments are available.

## ----eval= F------------------------------------------------------------------
# res <- ANOVA_DDCt(
#   data_3factor,
#   numOfFactors = 3,
#   numberOfrefGenes = 1,
#   specs = "Conc",
#   block = NULL)
# 
# model <- res$perGene$E_PO$lm
# # Relative expression values for Concentration main effect
# Means_DDCt(model, specs = "Conc")
# 
# # contrast        RE        SE df       LCL       UCL p.value sig
# # L vs H   0.1703610 0.2208988 24 0.1242014 0.2336757 <0.0001 ***
# # M vs H   0.2227247 0.2208988 24 0.1623772 0.3055004 <0.0001 ***
# # M vs L   1.3073692 0.2208988 24 0.9531359 1.7932535  0.0928 .
# #
# #Results are averaged over the levels of: Type, SA
# #Confidence level used: 0.95
# 
# # Relative expression values for Concentration sliced by Type
# Means_DDCt(model, specs = "Conc | Type")
# 
# #Type = R:
# # contrast       RE        SE df       LCL      UCL p.value sig
# # L vs H   0.103187 0.3123981 24 0.0659984 0.161331 <0.0001 ***
# # M vs H   0.339151 0.3123981 24 0.2169210 0.530255 <0.0001 ***
# # M vs L   3.286761 0.3123981 24 2.1022126 5.138776 <0.0001 ***
# #
# #Type = S:
# # contrast       RE        SE df       LCL      UCL p.value sig
# # L vs H   0.281265 0.3123981 24 0.1798969 0.439751 <0.0001 ***
# # M vs H   0.146266 0.3123981 24 0.0935518 0.228684 <0.0001 ***
# # M vs L   0.520030 0.3123981 24 0.3326112 0.813055  0.0059 **
# #
# #Results are averaged over the levels of: SA
# #Confidence level used: 0.95
# 
# # Relative expression values for Concentration sliced by Type and SA
# Means_DDCt(model, specs = "Conc | Type * SA")

## ----eval= F------------------------------------------------------------------
# data <- read.csv(system.file("extdata", "data_repeated_measure_1.csv", package = "rtpcr"))
# res3 <- ANOVA_DDCt(
#   data,
#   numOfFactors = 1,
#   numberOfrefGenes = 1,
#   specs = "time",
#   block = NULL,
#   model = wDCt ~ time + (1 | id))
# 
# residuals <- resid(res3$perGene$Target$lm)
# shapiro.test(residuals)
# par(mfrow = c(1,2))
# plot(residuals)
# qqnorm(residuals)
# qqline(residuals, col = "red")

## ----eval= F------------------------------------------------------------------
# # Example input data frame with technical replicates
# data1 <- read.csv(system.file("extdata", "data_withTechRep.csv", package = "rtpcr"))
# 
# # Calculate mean of technical replicates using first four columns as groups
# meanTech(data1,
#          groups = 1:2,
#          numOfFactors = 1,
#          block = NULL)

