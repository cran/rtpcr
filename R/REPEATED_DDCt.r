#' @title Fold change (\eqn{\Delta \Delta C_T} method) analysis of repeated measure qPCR data
#' 
#' @description \code{REPEATED_DDCt} function performs fold change (\eqn{\Delta \Delta C_T} method) 
#' analysis of observations repeatedly taken over different time courses. 
#' Data may be obtained over time from a uni- or multi-factorial experiment. The bar plot of the fold changes (FC) 
#' values along with the standard error (se) or confidence interval (ci) is also returned by the \code{REPEATED_DDCt} function. 
#' 
#' @details The \code{REPEATED_DDCt} function performs fold change (FC) analysis of observations repeatedly taken over time. 
#' The intended factor (could be time or any other factor) is defined for the analysis by the \code{factor} argument, 
#' then the function performs FC analysis on its levels
#' so that the first levels (as appears in the input data frame) is used as reference or calibrator. 
#' The function returns FC values along with confidence interval and standard error for the FC values.
#' 
#' @author Ghader Mirzaghaderi
#' @export REPEATED_DDCt
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import emmeans
#' @import lmerTest
#' @param x input data frame in which the first column is id, 
#' followed by the factor column(s) which include at least time. 
#' The first level of time in data frame is used as calibrator or reference level.
#' Additional factor(s) may also be present. Other columns are efficiency and Ct values of target and reference genes.
#'  \strong{NOTE:} In the "id" column, a unique number is assigned to each individual from which samples have been taken over time, 
#' for example see \code{data_repeated_measure_1}, 
#' all the three number 1 indicate one individual which has been sampled over three different time courses.
#' See \href{../doc/vignette.html}{\code{vignette}}, section "data structure and column arrangement" for details.
#' @param numberOfrefGenes number of reference genes which is 1 or 2 (Up to two reference genes can be handled).
#' as reference or calibrator which is the reference level or sample that all others are compared to. Examples are untreated 
#' of time 0. The FC value of the reference or calibrator level is 1 because it is not changed compared to itself.
#' If NULL, the first level of the main factor column is used as calibrator.
#' @param factor the factor for which the FC values is analysed. The first level of the specified factor in the input data frame is used as calibrator.
#' @param x.axis.labels.rename a vector replacing the x axis labels in the bar plot
#' @param block column name of the block if there is a blocking factor (for correct column arrangement see example data.). 
#' Block effect is usually considered as random and its interaction with any main effect is not considered.
#' @param p.adj Method for adjusting p values
#' @param plot  if \code{FALSE}, prevents the plot.
#' @param plotType  Plot based on "RE" (relative expression) or "log2FC" (log2 fold change).
#' @return A list with 5 elements:
#' \describe{
#'   \item{Final_data}{Input data frame plus the weighted Delat Ct values (wDCt)}
#'   \item{lm}{lm of factorial analysis-tyle}
#'   \item{ANOVA_table}{ANOVA table}
#'   \item{FC Table}{Table of FC values, significance, confidence interval and standard error with the lower and upper limits for the selected factor levels.}
#'   \item{Bar plot of FC values}{Bar plot of the fold change values for the main factor levels.}
#' }
#' 
#' 
#' 
#' @examples
#' 
#' REPEATED_DDCt(data_repeated_measure_1,
#'             numberOfrefGenes = 1,
#'             factor = "time", block = NULL)
#'
#' REPEATED_DDCt(data_repeated_measure_2,
#'              numberOfrefGenes = 1,
#'              factor = "time", block = NULL)
#'                                                        
#'                                                        



REPEATED_DDCt <- function(x, 
                          numberOfrefGenes,
                          factor, 
                          block,
                          x.axis.labels.rename = "none",
                          p.adj = "none",
                          plot = TRUE,
                          plotType = "RE"){
  
  
  if (missing(numberOfrefGenes)) {
    stop("argument 'numberOfrefGenes' is missing, with no default")
  }
  if (missing(factor)) {
    stop("argument 'factor' is missing, with no default")
  }
  if (missing(block)) {
    stop("argument 'block' is missing, with no default. Requires NULL or a blocking factor column.")
  }
  
  
  if (is.null(block)) {
    
    
    if(numberOfrefGenes == 1) {
      id <- colnames(x)[1]
      if((ncol(x)-5) <= 1){
        factors = NULL
      } else {
        factors <- colnames(x)[2:(ncol(x)-5)]
      }
      colnames(x)[ncol(x)-4] <- "time"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-(log2(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      id <- colnames(x)[1]
      factors <- colnames(x)[2:(ncol(x)-7)]
      colnames(x)[ncol(x)-6] <- "time"
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
      id <- colnames(x)[1]
      factors <- colnames(x)[2:(ncol(x)-6)]
      colnames(x)[ncol(x)-5] <- "block"
      colnames(x)[ncol(x)-4] <- "time"
      colnames(x)[ncol(x)-3] <- "Etarget"
      colnames(x)[ncol(x)-2] <- "Cttarget"
      colnames(x)[ncol(x)-1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget)*x$Cttarget)-(log2(x$Eref)*x$Ctref))
      
    } else if(numberOfrefGenes == 2) {
      id <- colnames(x)[1]
      factors <- colnames(x)[2:(ncol(x)-8)]
      colnames(x)[ncol(x)-7] <- "block"
      colnames(x)[ncol(x)-6] <- "time"
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
  
  for (i in 2:which(names(x) == "time")) {
    x[[i]] <- factor(x[[i]], levels = unique(x[[i]]))
  }
  
  
  # Check if there is block
  if (is.null(block)) {
    if((ncol(x)-5) <= 2){
      formula <- wDCt ~ time + (1 | id)
    } else {
      formula <- paste("wDCt ~", paste("time"," *"), paste(factors, collapse = " * "), "+ (1 | id)")
    }
  } else {
    if((ncol(x)-6) <= 2){
      formula <- wDCt ~ time + (1|id) + (1|block/id)
    } else {
      formula <- paste("wDCt ~ ", paste("time"," *"), paste(factors, collapse = " * "), "+ (1 | id) + (1|block/id)")
    }
  }
  lm <- lmerTest::lmer(formula, data = x)
  ANOVA <- stats::anova(lm) 
  
  
  
  
  
  v <- match(colnames(x), factor)
  n <- which(!is.na(v))
  factor <- colnames(x)[n]
  lvls <- unique(x[,n])
  calibrartor <- lvls[1]
  
  warning(paste("The level", calibrartor, " of the selected factor was used as calibrator."))
  pp1 <- emmeans(lm, factor, data = x, adjust = p.adj, mode = "satterthwaite")
  pp2 <- as.data.frame(graphics::pairs(pp1), adjust = p.adj)
  if (length(lvls) >= 3){
    pp3 <- pp2[1:length(lvls) - 1,] 
  } else {
    pp3 <- pp2
  }
  ci <- as.data.frame(stats::confint(graphics::pairs(pp1)), adjust = p.adj)[1:length(lvls)-1,]
  pp <- cbind(pp3, lower.CL = ci$lower.CL, upper.CL = ci$upper.CL)
  
  bwDCt <- x$wDCt   
  se <- summarise(
    group_by(data.frame(factor = x[n], bwDCt = bwDCt), x[n]),
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
  
  
  words <- strsplit(as.character(contrast[1]), " ")[[1]]
  referencelevel <- words[1]
 
  
  reference <- data.frame(contrast = as.character(referencelevel),
                          RE = 1,
                          log2FC = 0,
                          pvalue = 1, 
                          sig = " ",
                          LCL = 0,
                          UCL = 0,
                          se = se$se[1])
  
  tableC  <- rbind(reference, post_hoc_test)
  
  #round tableC to 4 decimal places
  tableC[, sapply(tableC, is.numeric)] <- lapply(tableC[, sapply(tableC, is.numeric)], function(x) round(x, 4))
  
  FINALDATA <- x
  
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
  FCp <- as.numeric(tableC$FC)
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
  

  outlist2 <- structure(list(Final_data = x,
                             lm = lm,
                             ANOVA_table = ANOVA,
                             FC_statistics_of_the_main_factor  = tableC,
                             RE_Plot_of_the_main_factor_levels = pfc1,
                             log2FC_Plot_of_the_main_factor_levels = pfc2), class = "XX")
  
  print.XX <- function(outlist2){
    print(outlist2$ANOVA_table)
    cat("\n", sep = '',"Fold Change table", "\n")
    print(outlist2$FC_statistics_of_the_main_factor)
    
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
