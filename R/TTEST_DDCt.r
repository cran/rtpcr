#' @title Expression (\eqn{\Delta \Delta C_T} method) analysis of target genes using t-test
#' 
#' @description t.test based analysis of the fold change expression for any number of target genes.
#' 
#' @details The \code{TTEST_DDCt} function applies a t.test based analysis to calculate 
#' fold change (\eqn{\Delta \Delta C_T} method) expression and returns related statistics for any number of 
#' target genes that have been evaluated under control and treatment conditions. This function also returns the 
#' expression bar plot based on fold change or log2 fold change. Sampling may be paired or unpaired. 
#' One or two reference genes can be used. Unpaired and paired samples are commonly analyzed 
#' using unpaired and paired t-test, respectively.  \strong{NOTE:} Paired samples in quantitative PCR refer to two sample 
#' data that are collected from one set of individuals 
#' at two different conditions, for example before and after a treatment or at two different time points. While 
#' for unpaired samples, two sets of individuals are used: one under untreated and the other set under treated 
#' condition.  Paired samples allow to compare gene expression changes within the same individual, reducing 
#' inter-individual variability. 
#' @author Ghader Mirzaghaderi
#' @export TTEST_DDCt
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @param x a data frame of 4 columns including Conditions, E (efficiency), Gene and Ct values (see examples below). Biological replicates needs to be equal for all Genes. Each Ct value is the mean of technical replicates. Complete amplification efficiencies of 2 is assumed here for all wells but the calculated efficienies can be used instead. See \href{../doc/vignette.html}{\code{vignette}} for details about "data structure and column arrangement".
#' @param paired a logical indicating whether you want a paired t-test.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param numberOfrefGenes number of reference genes. Up to two reference genes can be handled.
#' @param p.adj Method ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none") for adjusting p values.  
#' @param order a vector determining genes order on the output graph.
#' @param plotType  Plot based on "RE" (relative expression) or "log2FC" (log2 fold change).
#' @return A list of two elements:
#' \describe{
#'   \item{Result}{Output table including the Fold Change values, lower and upper confidence interval, pvalue and standard error with the lower and upper limits.}
#' }
#' For more information about the test procedure and its arguments,
#' refer \code{\link[stats]{t.test}}, and \code{\link[stats]{lm}}.
#' If the residuals of the model do not follow normal distribution and variances between the two groups are not homoGene, \code{\link[stats]{wilcox.test}} procedure may be concidered
#'
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method. Methods 25 (4). doi:10.1006/meth.2001.1262.
#'
#' Ganger, MT, Dietz GD, and Ewing SJ. 2017. A common base method for analysis
#' of qPCR data and the application of simple blocking in qPCR experiments.
#' BMC bioinformatics 18, 1-11.
#'
#' Yuan, Joshua S, Ann Reed, Feng Chen, and Neal Stewart. 2006.
#' Statistical Analysis of Real-Time PCR Data. BMC Bioinformatics 7 (85).
#' doi:10.1186/1471-2105-7-85.
#'
#' @examples
#'
#' # See the sample data structure
#' data_ttest
#'
#' # Getting t.test results
#' TTEST_DDCt(data_ttest,
#'    paired = FALSE,
#'    var.equal = TRUE,
#'    numberOfrefGenes = 1)
#'
#'
#'
#' TTEST_DDCt(Taylor_etal2019, 
#'           numberOfrefGenes = 2, 
#'           var.equal = TRUE,
#'           p.adj = "BH")
#'  
#'
#'




TTEST_DDCt <- function(x,
                       numberOfrefGenes, 
                       paired = FALSE, 
                       var.equal = TRUE, 
                       p.adj = "BH",
                       order = "none", 
                       plotType = "RE") {
  
  colnames(x)[1] <- "Condition"
  colnames(x)[2] <- "Gene"
  colnames(x)[3] <- "E"
  colnames(x)[4] <- "Ct"
  
  default.order <- unique(x[,2])[-length(unique(x[,2]))] 
  
  
  r <- nrow(x)/(2 * length(unique(x$Gene)))
  
  if(!all(r %% 1 == 0)) {
    stop("Error: Replicates are not equal for all Genes!")
  } else {
    
    
    
    x <- data.frame(x, wCt = log2(x$E) * x$Ct)
    
    if (numberOfrefGenes > 2) stop("Only up to 2 reference genes can be handled!")
    if(numberOfrefGenes == 1) {
      x <- x
    } else if (numberOfrefGenes == 2) {
      a <- (((2 * r) * (length(unique(x$Gene)) - 2)) + 1)
      b <- ((length(unique(x$Gene)) - 1) * 2 * r)
      mwCT <- (x$wCt[a:b] + x$wCt[(a+(2*r)):(b+(2*r))])/2
      x$wCt[a:b] <- mwCT
      x <- x[-((a+(2*r)):(b+(2*r))),]
    }
    
    
    
    GENE <- x$Gene
    
    
    levels_to_compare <- unique(GENE)[-length(unique(GENE))]
    res <- matrix(nrow = length(levels_to_compare), ncol = 7)
    colnames(res) <- c("Gene", "dif", "RE", "LCL", "UCL", "pvalue", "se")
    subset_df <- data.frame(group = character(), Gene = character(), wDCt = numeric())
    for (i in 1:length(levels_to_compare)) {
      subset_df <- rbind(subset_df, 
                         data.frame(group = x[GENE == levels_to_compare[i], "Condition"],
                                    Gene = levels_to_compare[i],
                                    wDCt = x[GENE == levels_to_compare[i], "wCt"] - x[GENE == utils::tail(unique(GENE), 1), "wCt"]))
    }
    
    
    
    
    subset <- matrix(NA, nrow = 2 * r, ncol=length(levels_to_compare))
    ttest_result <- vector("list", length(levels_to_compare))
    
    for (i in 1:length(levels_to_compare)) {
      subset[,i] <- x[GENE == levels_to_compare[i], "wCt"] - x[GENE == utils::tail(unique(GENE), 1), "wCt"]
      ttest_result[[i]] <- stats::t.test(subset[(r + 1):(2 * r), i], subset[1:r, i], paired = paired, var.equal = var.equal)
      
      res[i, ] <- c(levels_to_compare[i],
                    round(mean(subset[(r+1):(2*r),i]) - mean(subset[1:r, i]), 4),
                    round(2^-((mean(subset[(r+1):(2*r), i]) - mean(subset[1:r,i]))), 4),
                    round(2^(-ttest_result[[i]]$conf.int[2]), 4), # Lower error bar point
                    round(2^(-ttest_result[[i]]$conf.int[1]), 4), # Upper error bar point
                    round(ttest_result[[i]]$p.value, 4),
                    round(stats::sd(subset[(r+1):(2*r),i])/sqrt(r), 4))
      
    }
    res <- as.data.frame(res)
    res$RE <- as.numeric(res$RE)
    res$se <- as.numeric(res$se)
    res$dif <- NULL
    res <- data.frame(res, 
                      log2FC = log2(res$RE),
                      Lower.se.RE = round(2^(log2(res$RE) - res$se), 4), 
                      Upper.se.RE = round(2^(log2(res$RE) + res$se), 4),
                      p.adj = stats::p.adjust(res$pvalue, method = p.adj))
    

a <- data.frame(res, d = 0, Lower.se.log2FC = 0, Upper.se.log2FC = 0)
res$Lower.se
    for (i in 1:length(res$RE)) {
      if (res$RE[i] < 1) {
        a$Lower.se.log2FC[i] <- (res$Upper.se.RE[i]*log2(res$RE[i]))/res$RE[i]
        a$Upper.se.log2FC[i] <- (res$Lower.se.RE[i]*log2(res$RE[i]))/res$RE[i]
        a$d[i] <- (res$Upper.se.RE[i]*log2(res$RE[i]))/res$RE[i] - 0.2
      } else {
        a$Lower.se.log2FC[i] <- (res$Lower.se.RE[i]*log2(res$RE[i]))/res$RE[i]
        a$Upper.se.log2FC[i] <- (res$Upper.se.RE[i]*log2(res$RE[i]))/res$RE[i]
        a$d[i] <- (res$Upper.se.RE[i]*log2(res$RE[i]))/res$RE[i] + 0.2
      }
    }

res <- data.frame(res, 
                  Lower.se.log2FC = a$Lower.se.log2FC,
                  Upper.se.log2FC = a$Upper.se.log2FC)


######################
# Getting barplot:

df2 <- res

# Convert the column to a factor with specified levels
if(any(order == "none")){
  df2$Gene <- factor(df2$Gene, levels = default.order)
} else {
  df2$Gene <- factor(df2$Gene, levels = order)
}


# Order the data frame based on the specified column
df2 <- df2[order(df2$Gene), ]

Gene <- df2$Gene
Gene <- factor(Gene, levels = Gene)
Fold_Change <- df2$RE
Lower.Er <- df2$LCL
Upper.Er <- df2$UCL
pvalue <- as.numeric(df2$p.adj)
label <- .convert_to_character(pvalue)
se <- df2$se
Upper.se.log2FC <- df2$Upper.se.log2FC
Lower.se.log2FC <- df2$Lower.se.log2FC


df2 <- data.frame(df2, d = 0)
for (i in 1:length(df2$RE)) {
  if (df2$RE[i] < 1) {
    df2$d[i] <- (df2$Upper.se.RE[i]*log2(df2$RE[i]))/df2$RE[i] - 0.2
  } else {
    df2$d[i] <- (df2$Upper.se.RE[i]*log2(df2$RE[i]))/df2$RE[i] + 0.2
  }
}



if(plotType == "RE"){
  p <- ggplot(df2, aes(Gene, as.numeric(Fold_Change))) + 
    geom_col() +
    geom_errorbar(aes(ymin = Lower.se.RE, ymax=Upper.se.RE), width=0.1) +
    geom_text(aes(label = label, x = Gene,
                  y = Upper.se.RE + 0.2)) +
    ylab("Relative expression (DDCt)")
}
if(plotType == "log2FC"){
  p <- ggplot(df2, aes(Gene, as.numeric(log2FC))) +
    geom_col() +
    geom_errorbar(aes(ymin = Upper.se.log2FC, ymax=Lower.se.log2FC), width=0.1) +
    geom_text(aes(label = label, x = Gene,
                  y = d)) +
    ylab("log2FC")
}
######################

    
    Raw_df <- melt(subset, value.name = "wDCt")[-1]
    res <- list(Result = res, plot = p)
    return(res)
  }
}
