#' @title ANOVA of RE values based on CRD 
#' @description Analysis of Variance of relative efficiency (RE) values based on a completely randomized design (CRD). Even there are more than a factor in the experiment, it is still possible to apply CRD analysis on the factor-level combinations as treatments. Analysis of variance based on factorial design or analysis of covariance can be performed using \code{qpcrANCOVA} function.  
#' @details The \code{qpcrANOVA} function performs analysis of variance (ANOVA) of relative efficiency (RE) values based on a completely randomized design (CRD). 
#' It is suitable when relative expression (RE) analysis between different treatment combinations 
#' (in a Uni- or multi-factorial experiment) is desired. If there are more than a factor in the experiment, 
#' it is still possible to apply CRD analysis on the factor-level combinations as treatments. 
#' For this, a column of treatment combinations is made first as a grouping factor Fold change analysis based 
#' on factorial design or analysis of covariance for the can be performed using \link{qpcrANCOVA}.
#' @author Ghader Mirzaghaderi
#' @export qpcrANOVA
#' @import dplyr
#' @import tidyr
#' @import reshape2
#' @import ggplot2
#' @import lme4
#' @import agricolae
#' @param x A data frame consisting of condition columns, target gene efficiency (E), target Gene Ct, reference gene efficiency and reference gene Ct values, respectively. Each Ct in the data frame is the mean of technical replicates. Complete amplification efficiencies of 2 was assumed in the example data for all wells but the calculated efficienies can be used instead.
#' @param numberOfrefGenes number of reference genes (1 or 2). Up to two reference genes can be handled.
#' @param block column name of the blocking factor (for correct column arrangement see example data.). When a qPCR experiment is done in multiple qPCR plates, variation resulting from the plates may interfere with the actual amount of gene expression. One solution is to conduct each plate as a complete randomized block so that at least one replicate of each treatment and control is present on a plate. Block effect is usually considered as random and its interaction with any main effect is not considered.
#' @param p.adj Method for adjusting p values (see p.adjust)
#' @return A list with 5 elements:
#' \describe{
#'   \item{Final_data}{The row data plus weighed delta Ct (wDCt) values.}
#'   \item{lm}{The output of linear model analysis including ANOVA tables based on factorial experiment and completely randomized design (CRD).}
#'   \item{ANOVA_factorial}{ANOVA table based on factorial arrangement}
#'   \item{ANOVA_CRD}{ANOVA table based on CRD}
#'   \item{Result}{The result table including treatments and factors, RE (Relative Expression), LCL, UCL, letter grouping and standard deviation of relative expression.}
#'   \item{Post_hoc_Test}{Post hoc test of FC (Fold Change), pvalue, significance and confidence interval (LCL, UCL).}
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
#' @examples
#'
#' # If the data include technical replicates, means of technical replicates
#' # should be calculated first using meanTech function.
#'
#' # Applying ANOVA analysis
#' qpcrANOVA(
#'      data_3factor,
#'      numberOfrefGenes = 1,
#'      p.adj = "none")
#'
#'
#' qpcrANOVA(
#'     data_2factorBlock,
#'     block = "Block",
#'     numberOfrefGenes = 1)
#'
#'



qpcrANOVA <- function(x,
                      numberOfrefGenes,
                      block = NULL,
                      p.adj = c("none","holm","hommel", 
                                "hochberg", "bonferroni", "BH", "BY", "fdr")){
  
  
  resultx <- .addwDCt(x)
  x<- resultx$x
  factors <- resultx$factors
  # Check if there is block
  if (is.null(block)) {
    # Concatenate the columns using paste0
    x$T <- do.call(paste, c(x[1:(ncol(x)-6)], sep = ":"))
    x
    lm <- stats::lm(wDCt ~ T, x)
    anovaCRD <- stats::anova(lm)
    
  } else {
    # Concatenate the columns using paste0
    x$T <- do.call(paste, c(x[1:(ncol(x)-7)], sep = ":"))
    x
    lm <- stats::lm(wDCt ~ block + T, x)
    anovaCRD <- stats::anova(lm)
  }
  
  
  
  # Preparing final result table including letter grouping of the means for T
  g <- LSD.test(lm, "T", group = T, console = F, alpha = 0.05, p.adj = p.adj)$groups
  g <- g[rev(rownames(g)),] #order the result the way you want
  g$groups <- .invOrder(as.character(g$groups))
  mean <- LSD.test(lm, "T", group = T, console = F, alpha = 0.05, p.adj = p.adj)$means
  
  
  # Comparing mean pairs that also returns CI for T
  # Preparing final result table including letter grouping of the means
  meanPP <- LSD.test(lm, "T", group = F, console = F, alpha = 0.05, p.adj = p.adj)
  meanPairs <- meanPP$comparison
  ROWS <- rownames(meanPairs)
  diffs <- meanPairs$difference
  pval <- meanPairs$pvalue
  signif <- meanPairs$signif.
  ucl <- meanPairs$UCL
  lcl <- meanPairs$LCL
  Post_hoc_Testing <- data.frame(row.names = ROWS,
                                 FC = round(10^(-diffs), 4),
                                 pvalue = pval,
                                 signif. = signif,
                                 LCL = round(10^(-ucl), 4),
                                 UCL = round(10^(-lcl), 4))
  
  
  RowNames <- rownames(mean)
  mean$RowNames <- RowNames
  mean <- separate(mean, RowNames, into = factors, sep = ":", remove = T)
  mean <- mean[order(rownames(mean)),]
  g <- g[order(rownames(g)),]
  
  bwDCt <- 10^(-x$wDCt)
  sdRow <- summarise(
    group_by(data.frame(T = x$T, bwDCt = bwDCt), T),
    sd = sd(bwDCt, na.rm = TRUE))
  sd <- sdRow[order(sdRow$T),]
  
  Results <- data.frame(mean[,(ncol(mean)-2):ncol(mean)],
                        RE = round(10^(-mean$wDCt), 4),
                        LCL = round(10^(-mean$LCL), 4),
                        UCL = round(10^(-mean$UCL), 4),
                        letters = g$groups,
                        std = round(sd$sd, 4))
  
  
  # removing additional columns!
  if(length(factors) == 1) {
    Results <- Results[, -(1:2)]
    
  } else if(length(factors) == 2) {
    Results <- Results[, -1]
    
  } else if(length(factors) == 3) {
    Results <- Results
  }
  
  
  
  xx <- x[, -(ncol(x))] # Removing the last column of T
  
  
  outlist <- list(Final_data = xx,
                  lmCRD = lm,
                  ANOVA_CRD = anovaCRD,
                  Result = Results,
                  Post_hoc_Test = Post_hoc_Testing)
  
  
  return(outlist)
}
