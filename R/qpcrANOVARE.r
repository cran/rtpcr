#' @title Relative expression (\eqn{\Delta C_T} method) analysis using ANOVA 
#' 
#' @description Analysis of variance of relative expression (\eqn{\Delta C_T} method) values for 
#' all factor level combinations in the experiment in which the expression level of a 
#' reference gene is used as normalizer. 
#' 
#' @details The \code{qpcrANOVARE} function performs analysis of variance (ANOVA) of relative 
#' expression (RE) values for all factor level combinations as treatments using the expression 
#' level of a reference gene is used as normalizer. To get a reliable result, the expression of 
#' the reference gene needs to be constant across all test samples and it expression should not 
#' be affected by the experimental treatment under study.
#' 
#' @author Ghader Mirzaghaderi
#' @export qpcrANOVARE
#' @import dplyr
#' @import tidyr
#' @import reshape2
#' @import lmerTest
#' @import multcomp
#' @import multcompView
#' @param x a data frame consisting of condition columns, target gene efficiency (E), target Gene Ct, reference 
#' gene efficiency and reference gene Ct values, respectively. Each Ct in the data frame is the mean of 
#' technical replicates. Complete amplification efficiencies of 2 was assumed in the example data for 
#' all wells but the calculated efficienies can be used instead.  \strong{NOTE:} Each line belongs to a separate 
#' individual reflecting a non-repeated measure experiment). See \href{../doc/vignette.html}{\code{vignette}}, 
#' section "data structure and column arrangement" for details.
#' 
#' @param numberOfrefGenes number of reference genes (1 or 2). Up to two reference genes can be handled.
#' @param block column name of the blocking factor (for correct column arrangement see example data.). 
#' When a qPCR experiment is done in multiple qPCR plates, variation resulting from the plates may 
#' interfere with the actual amount of gene expression. One solution is to conduct each plate as a 
#' complete randomized block so that at least one replicate of each treatment and control is present 
#' on a plate. Block effect is usually considered as random and its interaction with any main effect is 
#' not considered.
#' @param alpha significance level
#' @param adjust method for adjusting p-values
#' @return A list with 4 elements:
#' \describe{
#'   \item{Final_data}{The row data plus weighed delta Ct (wDCt) values.}
#'   \item{lm}{The output of linear model analysis including ANOVA tables}
#'   \item{ANOVA}{ANOVA table based on CRD}
#'   \item{Result}{The result table including treatments and factors, RE (Relative Expression), LCL, UCL, 
#'   letter display for pair-wise comparisons and standard error with the lower and upper limits.}
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
#' # Applying ANOVA
#' qpcrANOVARE(data_3factor, numberOfrefGenes = 1, block = NULL)
#'
#'
#' qpcrANOVARE(data_2factorBlock, block = "Block", numberOfrefGenes = 1)
#'
#'


qpcrANOVARE <- function(x, numberOfrefGenes, block, alpha = 0.05, adjust= "none")
{
  
  
  if (missing(numberOfrefGenes)) {
    stop("argument 'numberOfrefGenes' is missing, with no default")
  }
  if (missing(block)) {
    stop("argument 'block' is missing, with no default. Requires NULL or a blocking factor column.")
  }
  
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
  
  
  # Check if there is block
  if(numberOfrefGenes == 1) {
    if (is.null(block)) {
      # Concatenate the columns using paste0
      x$T <- do.call(paste, c(x[1:(ncol(x)-6)], sep = ":"))
      x$T <- as.factor(x$T)
      lm <- stats::lm(wDCt ~ T, x)
      anovaCRD <- stats::anova(lm)
      
    } else {
      # Concatenate the columns using paste0
      x$T <- do.call(paste, c(x[1:(ncol(x)-7)], sep = ":"))
      x$T <- as.factor(x$T)
      lm <- stats::lm(wDCt ~ block + T, x)
      anovaCRD <- stats::anova(lm)
    }
  } 
  if(numberOfrefGenes == 2) {
    if (is.null(block)) {
      # Concatenate the columns using paste0
      x$T <- do.call(paste, c(x[1:(ncol(x)-8)], sep = ":"))
      x$T <- as.factor(x$T)
      lm <- stats::lm(wDCt ~ T, x)
      anovaCRD <- stats::anova(lm)
      
    } else {
      # Concatenate the columns using paste0
      x$T <- do.call(paste, c(x[1:(ncol(x)-9)], sep = ":"))
      x$T <- as.factor(x$T)
      lm <- stats::lm(wDCt ~ block + T, x)
      anovaCRD <- stats::anova(lm)
    }
  }
  
  
  emg = emmeans(lm,pairwise ~ T)
  meanPairs <- cld(emg[[1]], adjust= adjust, alpha = alpha, reversed = FALSE, Letters = letters)
  ROWS <- meanPairs$T
  diffs <- meanPairs$emmean
  ucl <- meanPairs$upper.CL
  lcl <- meanPairs$lower.CL
  letters <- meanPairs$.group 
  
  bwDCt <- x$wDCt 
  sdRow <- summarise(
    group_by(data.frame(T = x$T, bwDCt = bwDCt), T),
    mean = mean(bwDCt, na.rm = TRUE), 
    se = stats::sd(bwDCt/sqrt(length(bwDCt)), na.rm = TRUE)) 
  se <- sdRow[order(sdRow$mean),] 
  
  Results <- data.frame(row.names = ROWS,
                        RE = round(2^(-diffs), 4),
                        LCL = round(2^(-ucl), 4),
                        UCL = round(2^(-lcl), 4),
                        se = round(se$se, 4),
                        letters)
  
  
  RowNames <- rownames(Results)
  Results$RowNames <- RowNames
  
  
  mean <- separate(Results, RowNames, into = factors, sep = ":", remove = T)
  rownames(mean) <- NULL
  Results <- mean %>% 
    select(-1:-5) %>%  # Select all columns except the first 5
    cbind(mean[, 1:5]) 
  
  
  xx <- x[, -(ncol(x))] # Removing the last column of T
  
  
  
  Lower.se <- round(2^(log2(Results$RE) - Results$se), 4) 
  Upper.se <- round(2^(log2(Results$RE) + Results$se), 4)
  Results <- data.frame(Results, 
                       Lower.se = Lower.se, 
                       Upper.se = Upper.se)
  
  Results <- Results %>%
    relocate(Lower.se, Upper.se, .before = letters)
  
  outlist2 <- structure(list(Final_data = xx,
                             lmCRD = lm,
                             ANOVA = anovaCRD,
                             Results = Results), class = "XX")
  
  print.XX <- function(outlist2){
    print(outlist2$ANOVA)
    cat("\n", sep = '',"Relative expression table", "\n")
    print(outlist2$Results)
    invisible(outlist2)
  }
  print.XX(outlist2)
}