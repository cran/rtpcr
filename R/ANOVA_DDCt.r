#' @title Relative expression (\eqn{\Delta \Delta Ct} method) analysis using ANOVA and ANCOVA
#' 
#' @description Relative expression analysis using \eqn{\Delta \Delta Ct} method can be done by 
#' ANOVA (analysis of variance) and ANCOVA (analysis of covariance) through the \code{ANOVA_DDCt} function, for uni- or multi-factorial experiment data. The bar plot of the relative expression (RE) 
#' values along with the standard error (se) and confidence interval (ci) is returned by 
#' the \code{ANOVA_DDCt} function. 
#' 
#' @details ANOVA (analysis of variance) and ANCOVA (analysis of covariance) analysis of relative expression 
#' using \eqn{\Delta \Delta Ct} method can be done using \code{ANOVA_DDCt} function, for uni- or multi-factorial experiment data. 
#' If there are more than one factor, RE value calculations for 
#' the `mainFactor.column` and the statistical analysis is performed based on a full model factorial 
#' experiment by default. However, if `ancova` is defined for the `analysisType` argument,
#' RE values are calculated for the levels of the `mainFactor.column` and the other factors are 
#' used as covariate(s) in the analysis of variance, but we should consider ANCOVA table:
#' if the interaction between the main factor and covariate is significant, ANCOVA is not appropriate in this case. 
#' ANCOVA is basically used when a factor is affected by uncontrolled quantitative covariate(s). 
#' For example, suppose that wDCt of a target gene in a plant is affected by temperature. The gene may 
#' also be affected by drought. Since we already know that temperature affects the target gene, we are 
#' interested to know if the gene expression is also altered by the drought levels. We can design an 
#' experiment to understand the gene behavior at both temperature and drought levels at the same time. 
#' The drought is another factor (the covariate) that may affect the expression of our gene under the 
#' levels of the first factor i.e. temperature. The data of such an experiment can be analyzed by ANCOVA 
#' or using ANOVA based on a factorial experiment using \code{ANOVA_DDCt}. \code{ANOVA_DDCt} function performs RE 
#' analysis even there is only one factor (without covariate or factor  variable). Bar plot of relative expression 
#' (RE) values along with the standard errors are also returned by the \code{ANOVA_DDCt} function.
#'  
#' @author Ghader Mirzaghaderi
#' @export ANOVA_DDCt
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import lmerTest
#' @import emmeans
#' @param x a data frame of condition(s), biological replicates, efficiency (E) and Ct values of 
#' target and reference genes. Each Ct value in the data frame is the mean of technical replicates. 
#' \strong{NOTE:} Each line belongs to a separate individual reflecting a non-repeated measure experiment). 
#' See \href{../doc/vignette.html}{\code{vignette}}, section "data structure and column arrangement" for details.
#' 
#' @param numberOfrefGenes number of reference genes which is 1 or 2 (up to two reference genes can be handled).
#' @param analysisType should be one of "ancova" or "anova". Default is "anova".
#' @param mainFactor.column the factor for which relative expression is calculated for its levels. 
#' If \code{ancova} is selected as \code{analysisType}, the remaining factors (if any) are considered as covariate(s).
#' @param mainFactor.level.order  NULL (default) or a vector of main factor level names. If \code{NULL}, 
#' the first level of the \code{mainFactor.column} is used 
#' as calibrator. If a vector of main factor levels (in any order) is specified, the first level in the vector is 
#' used as calibrator. Calibrator is the reference level or sample that all others are compared to. Examples are untreated 
#' of time 0. The RE value of the reference or calibrator level is 1 because it is not changed compared to itself.
#' 
#' @param x.axis.labels.rename a vector replacing the x axis labels in the bar plot
#' @param block column name of the block if there is a blocking factor (for correct column arrangement see 
#' example data.). When a qPCR experiment is done in multiple qPCR plates, variation resulting from the 
#' plates may interfere with the actual amount of gene expression. One solution is to conduct each plate 
#' as a complete randomized block so that at least one replicate of each treatment and control is present 
#' on a plate. Block effect is usually considered as random and its interaction with any main effect 
#' is not considered.
#' @param p.adj Method for adjusting p values
#' @param plot  if \code{FALSE}, prevents the plot.
#' @param plotType  Plot based on "RE" (relative expression) or "log2FC" (log2 fold change).
#' @return A list with 7 elements:
#' \describe{
#'   \item{Final_data}{Input data frame plus the weighted Delat Ct values (wDCt)}
#'   \item{lm_ANOVA}{lm of factorial analysis-tyle}
#'   \item{lm_ANCOVA}{lm of ANCOVA analysis-type}
#'   \item{ANOVA_table}{ANOVA table}
#'   \item{ANCOVA_table}{ANCOVA table}
#'   \item{RE Table}{Table of RE (relative expression) values, log2FC (log2 fold change) values, significance and confidence interval and standard error with the lower and upper limits for the main factor levels.}
#'   \item{Bar plot of RE values}{Bar plot of the relative expression values for the main factor levels.}
#' }
#' 
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method. Methods 25 (4). doi:10.1006/meth.2001.1262.
#'
#' Ganger, MT, Dietz GD, and Ewing SJ. 2017. A common base method for analysis of qPCR data
#' and the application of simple blocking in qPCR experiments. BMC bioinformatics 18, 1-11.
#'
#' Yuan, Joshua S, Ann Reed, Feng Chen, and Neal Stewart. 2006.
#' Statistical Analysis of Real-Time PCR Data. BMC Bioinformatics 7 (85). doi:10.1186/1471-2105-7-85.
#' 
#' 
#' 
#' @examples
#'
#'
#' ANOVA_DDCt(data_1factor, numberOfrefGenes = 1, mainFactor.column = 1, block = NULL)
#' 
#'
#' ANOVA_DDCt(data_2factor, 
#' numberOfrefGenes = 1, 
#' mainFactor.column = 2, block = NULL, 
#' analysisType = "ancova")
#'
#'
#' # Data from Lee et al., 2020, Here, the data set contains technical 
#' # replicates so we get mean of technical replicates first:
#' df <- meanTech(Lee_etal2020qPCR, groups = 1:3)
#' ANOVA_DDCt(df, numberOfrefGenes = 1, analysisType = "ancova", block = NULL, 
#' mainFactor.column = 2, plotType = "log2FC")
#' 
#' 
#' ANOVA_DDCt(data_2factorBlock,  
#'     numberOfrefGenes = 1, 
#'     mainFactor.column = 1, 
#'     mainFactor.level.order = c("S", "R"), 
#'     block = "block", 
#'     analysisType = "ancova")
#'
#' 
#' 
#' df <- meanTech(Lee_etal2020qPCR, groups = 1:3) 
#' df2 <- df[df$factor1 == "DSWHi",][-1]
#' ANOVA_DDCt(df2,
#' mainFactor.column = 1,
#' block = NULL, 
#' numberOfrefGenes = 1, 
#' analysisType = "anova")
#' 
#'
#' addline_format <- function(x,...){gsub('\\s','\n',x)}
#' ANOVA_DDCt(data_1factor, numberOfrefGenes = 1, 
#' mainFactor.column = 1,
#' block = NULL, 
#' x.axis.labels.rename = addline_format(c("Control", 
#'                                         "Treatment_1 vs Control", 
#'                                         "Treatment_2 vs Control")))
#'                                                        
#'                                                        



ANOVA_DDCt <- function(
    x, 
    numberOfrefGenes = 1, 
    mainFactor.column = 1, 
    analysisType = "anova",
    mainFactor.level.order = NULL, 
    block = NULL, 
    x.axis.labels.rename = "none",
    p.adj = "none",  
    plot = TRUE,
    plotType = "RE"
)
{
  
  
  
  if (missing(numberOfrefGenes)) {
    stop("argument 'numberOfrefGenes' is missing, with no default")
  }
  if (missing(mainFactor.column)) {
    stop("argument 'mainFactor.column' is missing, with no default")
  }
  if (missing(block)) {
    stop("argument 'block' is missing, with no default. Requires NULL or a blocking factor column.")
  }
  
  
  x <- x[, c(mainFactor.column, (1:ncol(x))[-mainFactor.column])] 
  
  
  if (is.null(mainFactor.level.order)) {
    x[,1] <- factor(x[,1], levels = unique(x[,1]))
    mainFactor.level.order <- unique(x[,1])
    calibrartor <- x[,1][1]
    #warning(paste("The", calibrartor, "level was used as calibrator."))
    warning(structure(paste("The", calibrartor, "level was used as calibrator."), foreground = "blue"))
  } else if (any(is.na(match(unique(x[,1]), mainFactor.level.order))) == TRUE){
    stop("The `mainFactor.level.order` doesn't match main factor levels.")
  } else {
    x <- x[order(match(x[,1], mainFactor.level.order)), ]
    x[,1] <- factor(x[,1], levels = mainFactor.level.order)
  }
  
  
  # The data frame doesn't have ? columns. Please refer to the vignette to ensure that our data is properly structured.
  
  
  if (is.null(block)) {
    
    
    if(numberOfrefGenes == 1) {
      
      factors <- colnames(x)[1:(ncol(x)-5)]
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-(log2(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      
      factors <- colnames(x)[1:(ncol(x)-7)]
      colnames(x)[ncol(x)-6] <- "rep"
      colnames(x)[ncol(x)-5] <- "Etarget"
      colnames(x)[ncol(x)-4] <- "Cttarget"
      colnames(x)[ncol(x)-3] <- "Eref"
      colnames(x)[ncol(x)-2] <- "Ctref"
      colnames(x)[ncol(x)-1] <- "Eref2"
      colnames(x)[ncol(x)] <- "Ctref2"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-
                        ((log2(x$Eref)*x$Ctref) + (log2(x$Eref2)*x$Ctref2))/2)
    }
    
  } else {
    if(numberOfrefGenes == 1) {
      
      factors <- colnames(x)[1:(ncol(x)-6)]
      colnames(x)[ncol(x)-5] <- "block"
      colnames(x)[ncol(x)-4] <- "rep"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-(log2(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      factors <- colnames(x)[1:(ncol(x)-8)]
      colnames(x)[ncol(x)-7] <- "block"
      colnames(x)[ncol(x)-6] <- "rep"
      colnames(x)[ncol(x)-5] <- "Etarget"
      colnames(x)[ncol(x)-4] <- "Cttarget"
      colnames(x)[ncol(x)-3] <- "Eref"
      colnames(x)[ncol(x)-2] <- "Ctref"
      colnames(x)[ncol(x)-1] <- "Eref2"
      colnames(x)[ncol(x)] <- "Ctref2"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-
                        ((log2(x$Eref)*x$Ctref) + (log2(x$Eref2)*x$Ctref2))/2)
    }
  }
  
  
  
  # converting columns 1 to time as factor
  
  for (i in 2:which(names(x) == "rep")-1) {
    x[[i]] <- factor(x[[i]], levels = unique(x[[i]]))
  }
  
  
  # Check if there is block
  if (is.null(block)) {
    
    # ANOVA based on factorial design
    formula_ANOVA <- paste("wDCt ~", paste(factors, collapse = " * "), "+ (1 | rep)")
    base::suppressMessages(lmf <- lmerTest::lmer(formula_ANOVA, data = x))
    ANOVA <- stats::anova(lmf)
    # ANCOVA 
    formula_ANCOVA <- paste("wDCt ~", paste(rev(factors), collapse = " + "), "+ (1 | rep)")
    base::suppressMessages(lmc <- lmerTest::lmer(formula_ANCOVA, data = x))
    ANCOVA <- stats::anova(lmc)
    
  } else {
    # If ANOVA based on factorial design was desired with blocking factor:
    formula_ANOVA <- paste("wDCt ~ block +", paste(factors, collapse = " * "), "+ (1 | rep)")
    lmfb <- lmerTest::lmer(formula_ANOVA, data = x)
    ANOVA <- stats::anova(lmfb)
    # ANCOVA 
    formula_ANCOVA <- paste("wDCt ~ block +", paste(rev(factors), collapse = " + "), "+ (1 | rep)")
    lmcb <- lmerTest::lmer(formula_ANCOVA, data = x)
    ANCOVA <- stats::anova(lmcb)
  }
  
  
  
  
  
  # Type of analysis: ancova or anova
  if (is.null(block)) {
    if(analysisType == "ancova") {
      lm <- lmc
    } 
    else{
      lm <- lmf
    }
  } else {
    if(analysisType == "ancova") {
      lm <- lmcb
    } 
    else{
      lm <- lmfb
    } 
  }
  
  
  
  
  
  
  
  
  
  pp1 <- emmeans(lm, colnames(x)[1], data = x, adjust = p.adj, mode = "satterthwaite")
  pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
  pp3 <- pp2[1:length(mainFactor.level.order)-1,]
  ci <- as.data.frame(stats::confint(graphics::pairs(pp1)), adjust = p.adj)[1:length(unique(x[,1]))-1,]
  pp <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)
  
  
  
  bwDCt <- x$wDCt   
  se <- summarise(
    group_by(data.frame(factor = x[,1], bwDCt = bwDCt), x[,1]),
    se = stats::sd(bwDCt, na.rm = TRUE)/sqrt(length(bwDCt)))  
  
  
  sig <- .convert_to_character(pp$p.value)
  contrast <- pp$contrast
  post_hoc_test <- data.frame(contrast, 
                              RE = 1/(2^-(pp$estimate)),
                              log2FC = log2(1/(2^-(pp$estimate))),
                              pvalue = pp$p.value,
                              sig = sig,
                              LCL = 1/(2^-pp$lower.CL),
                              UCL = 1/(2^-pp$upper.CL),
                              se = se$se[-1])
  
  reference <- data.frame(contrast = mainFactor.level.order[1],
                          RE = 1,
                          log2FC = 0,
                          pvalue = 1, 
                          sig = " ",
                          LCL = 0,
                          UCL = 0,
                          se = se$se[1])
  
  tableC <- rbind(reference, post_hoc_test)
  
  #round tableC to 4 decimal places
  tableC[, sapply(tableC, is.numeric)] <- lapply(tableC[, sapply(tableC, is.numeric)], function(x) round(x, 4))
  
  FINALDATA <- x
  
  tableC$contrast <- as.character(tableC$contrast)
  tableC$contrast <- sapply(strsplit(tableC$contrast, " - "), function(x) paste(rev(x), collapse = " vs "))
  
  if(any(x.axis.labels.rename == "none")){
    tableC
  }else{
    tableC$contrast <- x.axis.labels.rename
  }
  
  
  
  
  tableC$contrast <- factor(tableC$contrast, levels = unique(tableC$contrast))
  contrast <- tableC$contrast
  LCL <- tableC$LCL
  UCL <- tableC$UCL
  REp <- as.numeric(tableC$RE)
  FCp <- as.numeric(tableC$log2FC)
  significance <- tableC$sig
  se <- tableC$se
  
  tableC <- data.frame(tableC, 
                       Lower.se.RE = round(2^(log2(tableC$RE) - tableC$se), 4), 
                       Upper.se.RE = round(2^(log2(tableC$RE) + tableC$se), 4))  
  ##################################################
  a <- data.frame(tableC, d = 0)

  for (i in 1:length(tableC$RE)) {
    if (tableC$RE[i] < 1) {
      a$Lower.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] - 0.2
    } else {
      a$Lower.se[i] <- (tableC$Lower.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$Upper.se[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i]
      a$d[i] <- (tableC$Upper.se.RE[i]*log2(tableC$RE[i]))/tableC$RE[i] + 0.2
    }
  }
  pfc1 <- ggplot(a, aes(contrast,RE)) + 
    geom_col() +
    geom_errorbar(aes(ymin = tableC$Lower.se.RE, ymax=tableC$Upper.se.RE), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = tableC$Upper.se.RE + 0.2)) +
    ylab("Relative Expression (DDCt)")
  pfc2 <- ggplot(a, aes(contrast,log2FC)) +
    geom_col() +
    geom_errorbar(aes(ymin = Upper.se, ymax=Lower.se), width=0.1) +
    geom_text(aes(label = sig, x = contrast,
                  y = d)) +
    ylab("log2FC")

  tableC <- data.frame(tableC, Lower.se.log2FC = a$Lower.se, Upper.se.log2FC = a$Upper.se)
  ##################################################  
  
  if (is.null(block)) {
    lm_ANOVA <- lmf
    lm_ANCOVA <- lmc
  } else {
    lm_ANOVA <- lmfb
    lm_ANCOVA <- lmcb
  }
  
  
  outlist2 <- structure(list(Final_data = x,
                             lm_ANOVA = lm_ANOVA,
                             lm_ANCOVA = lm_ANCOVA,
                             ANOVA_table = ANOVA,
                             ANCOVA_table = ANCOVA,
                             Fold_Change  = tableC,
                             RE_Plot_of_the_main_factor_levels = pfc1,
                             log2FC_Plot_of_the_main_factor_levels = pfc2), class = "XX")
  
  print.XX <- function(outlist2){
    cat("ANOVA table", "\n")
    print(outlist2$ANOVA_table)
    cat("\n", sep = '', "ANCOVA table", "\n")
    print(outlist2$ANCOVA_table)
    cat("\n", sep = '', "Expression table", "\n")
    print(outlist2$Fold_Change)
    
    
    if (plot == TRUE){
      if(plotType == "RE"){
        cat("\n", sep = '', "Expression plot of main factor levels", "\n")
        print(outlist2$RE_Plot_of_the_main_factor_levels)
      }else{
        cat("\n", sep = '', "Expression plot of main factor levels", "\n")
        print(outlist2$log2FC_Plot_of_the_main_factor_levels)
      }
    }
    invisible(outlist2)
  }
  print.XX(outlist2)
}
