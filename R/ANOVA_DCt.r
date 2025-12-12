#' @title Relative expression (\eqn{\Delta Ct} method) analysis using ANOVA 
#' 
#' @description Analysis of variance of relative expression (\eqn{\Delta Ct} method) values for 
#' all factor level combinations in which the expression level of a 
#' reference gene is used as normalizer. 
#' 
#' @details The \code{ANOVA_DCt} function performs analysis of variance (ANOVA) of relative 
#' expression (RE) values for all factor level combinations using the expression 
#' level of reference gene(s) as a normalizer.
#' 
#' @author Ghader Mirzaghaderi
#' @export ANOVA_DCt
#' @import dplyr
#' @import tidyr
#' @import reshape2
#' @import lmerTest
#' @import multcomp
#' @import emmeans
#' @import multcompView
#' @param x a data.frame structured as described in the vignette consisting of condition columns, target gene efficiency (E), target Gene Ct, reference 
#' gene efficiency and reference gene Ct values, respectively. Each Ct in the data frame is the mean of 
#' technical replicates. Complete amplification efficiencies of 2 was assumed in the example data for 
#' all wells but the calculated efficienies can be used instead.  \strong{NOTE:} Each line belongs to a separate 
#' individual reflecting a non-repeated measure experiment). See \href{../doc/vignette.html}{\code{vignette}}, 
#' section "data structure and column arrangement" for details.
#' 
#' @param numberOfrefGenes number of reference genes (1 or 2). Up to two reference genes can be handled.
#' @param block A string or NULL. If provided, this should be the name of the column in \code{x} that 
#' indicates the blocking factor. When qPCR is performed on different plates, variations between plates can introduce noise 
#' that obscures the true biological differences in gene expression. By using a blocking factor, 
#' each treatment and control group can be replicated across different plates, ensuring that 
#' at least one replicate of each condition is present on every plate. 
#' @param alpha significance level for cld (default 0.05)
#' @param adjust p-value adjustment method passed to emmeans/cld
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
#' ANOVA_DCt(data_3factor, numberOfrefGenes = 1, block = NULL)
#'
#'
#' ANOVA_DCt(data_2factorBlock, block = "Block", numberOfrefGenes = 1)
#'
#'


ANOVA_DCt <- function(x,
                      numberOfrefGenes,
                      block,
                      alpha = 0.05,
                      adjust = "none") {
  
  ## ---- basic argument checks ----
  if (missing(numberOfrefGenes)) {
    stop("argument 'numberOfrefGenes' is missing, with no default")
  }
  if (missing(block)) {
    stop("argument 'block' is missing, with no default. Requires NULL or a blocking factor column.")
  }
  
  ## ---- rename columns and compute wDCt ----
  if (is.null(block)) {
    
    if (numberOfrefGenes == 1) {
      factors <- colnames(x)[1:(ncol(x) - 5)]
      colnames(x)[ncol(x) - 4] <- "rep"
      colnames(x)[ncol(x) - 3] <- "Etarget"
      colnames(x)[ncol(x) - 2] <- "Cttarget"
      colnames(x)[ncol(x) - 1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget) * x$Cttarget) - (log2(x$Eref) * x$Ctref))
      
    } else if (numberOfrefGenes == 2) {
      factors <- colnames(x)[1:(ncol(x) - 7)]
      colnames(x)[ncol(x) - 6] <- "rep"
      colnames(x)[ncol(x) - 5] <- "Etarget"
      colnames(x)[ncol(x) - 4] <- "Cttarget"
      colnames(x)[ncol(x) - 3] <- "Eref"
      colnames(x)[ncol(x) - 2] <- "Ctref"
      colnames(x)[ncol(x) - 1] <- "Eref2"
      colnames(x)[ncol(x)] <- "Ctref2"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget) * x$Cttarget) -
                        ((log2(x$Eref) * x$Ctref) + (log2(x$Eref2) * x$Ctref2)) / 2)
    }
    
  } else { # with block
    
    if (numberOfrefGenes == 1) {
      factors <- colnames(x)[1:(ncol(x) - 6)]
      colnames(x)[ncol(x) - 5] <- "block"
      colnames(x)[ncol(x) - 4] <- "rep"
      colnames(x)[ncol(x) - 3] <- "Etarget"
      colnames(x)[ncol(x) - 2] <- "Cttarget"
      colnames(x)[ncol(x) - 1] <- "Eref"
      colnames(x)[ncol(x)] <- "Ctref"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget) * x$Cttarget) - (log2(x$Eref) * x$Ctref))
      
    } else if (numberOfrefGenes == 2) {
      factors <- colnames(x)[1:(ncol(x) - 8)]
      colnames(x)[ncol(x) - 7] <- "block"
      colnames(x)[ncol(x) - 6] <- "rep"
      colnames(x)[ncol(x) - 5] <- "Etarget"
      colnames(x)[ncol(x) - 4] <- "Cttarget"
      colnames(x)[ncol(x) - 3] <- "Eref"
      colnames(x)[ncol(x) - 2] <- "Ctref"
      colnames(x)[ncol(x) - 1] <- "Eref2"
      colnames(x)[ncol(x)] <- "Ctref2"
      
      x <- data.frame(x, wDCt = (log2(x$Etarget) * x$Cttarget) -
                        ((log2(x$Eref) * x$Ctref) + (log2(x$Eref2) * x$Ctref2)) / 2)
    }
  }
  
  ## ---- build treatment factor T and fit lm ----
  if (numberOfrefGenes == 1) {
    if (is.null(block)) {
      x$T <- do.call(paste, c(x[1:(ncol(x) - 6)], sep = ":"))
      x$T <- as.factor(x$T)
      lm_fit <- stats::lm(wDCt ~ T, data = x)
      anovaCRD <- stats::anova(lm_fit)
    } else {
      x$T <- do.call(paste, c(x[1:(ncol(x) - 7)], sep = ":"))
      x$T <- as.factor(x$T)
      lm_fit <- stats::lm(wDCt ~ block + T, data = x)
      anovaCRD <- stats::anova(lm_fit)
    }
  } else if (numberOfrefGenes == 2) {
    if (is.null(block)) {
      x$T <- do.call(paste, c(x[1:(ncol(x) - 8)], sep = ":"))
      x$T <- as.factor(x$T)
      lm_fit <- stats::lm(wDCt ~ T, data = x)
      anovaCRD <- stats::anova(lm_fit)
    } else {
      x$T <- do.call(paste, c(x[1:(ncol(x) - 9)], sep = ":"))
      x$T <- as.factor(x$T)
      lm_fit <- stats::lm(wDCt ~ block + T, data = x)
      anovaCRD <- stats::anova(lm_fit)
    }
  } else {
    stop("numberOfrefGenes must be 1 or 2")
  }
  
  ## ---- emmeans / multiple comparisons ----
  emg <- emmeans::emmeans(lm_fit, pairwise ~ T, mode = "satterthwaite")
  # use cld() on the emmeans object (the first element)
  meanPairs <- multcomp::cld(emg[[1]], adjust = adjust, alpha = alpha, reversed = FALSE, Letters = letters)
  # meanPairs typically contains columns: T, emmean, lower.CL, upper.CL, .group
  ROWS <- as.character(meanPairs[[1]])     # treatment labels in the same order
  diffs <- meanPairs$emmean
  ucl <- meanPairs$upper.CL
  lcl <- meanPairs$lower.CL
  letters_grp <- meanPairs$.group
  
  ## ---- compute group-wise means and SE (base R, robust) ----
  bwDCt <- x$wDCt
  means_by_T <- tapply(bwDCt, x$T, function(z) mean(z, na.rm = TRUE))
  sds_by_T   <- tapply(bwDCt, x$T, function(z) stats::sd(z, na.rm = TRUE))
  n_by_T     <- tapply(bwDCt, x$T, function(z) sum(!is.na(z)))
  se_by_T    <- sds_by_T / sqrt(n_by_T)
  
  se_df <- data.frame(T = names(means_by_T),
                      mean = as.numeric(means_by_T),
                      se = as.numeric(se_by_T),
                      stringsAsFactors = FALSE)
  
  # match se to the order used by emmeans/cld (ROWS)
  se_matched <- se_df$se[match(ROWS, se_df$T)]
  
  ## ---- build Results table ----
  Results <- data.frame(row.names = ROWS,
                        RE = round(2^(-diffs), 4),
                        log2FC = log2(round(2^(-diffs), 4)),
                        LCL = round(2^(-lcl), 4),
                        UCL = round(2^(-ucl), 4),
                        se = round(se_matched, 4),
                        letters = trimws(letters_grp), 
                        stringsAsFactors = FALSE)
  
  # preserve rownames as a column for splitting
  Results$RowNames <- rownames(Results)


  ## ---- split RowNames back to factor columns (base R) ----
  parts <- strsplit(Results$RowNames, ":", fixed = TRUE)
  parts_mat <- do.call(rbind, lapply(parts, function(p) {
    # ensure length matches number of factor columns
    length(p) <- length(factors)
    p
  }))
  parts_df <- as.data.frame(parts_mat, stringsAsFactors = FALSE)
  names(parts_df) <- factors
  
  # combine parts_df (factor columns) with Results
  Results_combined <- cbind(parts_df, Results)
  rownames(Results_combined) <- NULL
  
  ## ---- compute Lower/Upper SE for RE and attach to Results ----
  Results_combined$Lower.se.RE <- round(2^(log2(Results_combined$RE) - Results_combined$se), 4)
  Results_combined$Upper.se.RE <- round(2^(log2(Results_combined$RE) + Results_combined$se), 4)
  
  ## ---- compute log2FC SE bounds (vectorized) ----
  # initialize
  Results_combined$Lower.se.log2FC <- 0
  Results_combined$Upper.se.log2FC <- 0
  
  # vectorized computation replacing the for loop
  idx_less1 <- Results_combined$RE < 1
  idx_ge1   <- !idx_less1
  
  # when RE < 1
  Results_combined$Lower.se.log2FC[idx_less1] <- (Results_combined$Upper.se.RE[idx_less1] *
                                                    log2(Results_combined$RE[idx_less1])) / Results_combined$RE[idx_less1]
  Results_combined$Upper.se.log2FC[idx_less1] <- (Results_combined$Lower.se.RE[idx_less1] *
                                                    log2(Results_combined$RE[idx_less1])) / Results_combined$RE[idx_less1]
  
  # when RE >= 1
  Results_combined$Lower.se.log2FC[idx_ge1] <- (Results_combined$Lower.se.RE[idx_ge1] *
                                                  log2(Results_combined$RE[idx_ge1])) / Results_combined$RE[idx_ge1]
  Results_combined$Upper.se.log2FC[idx_ge1] <- (Results_combined$Upper.se.RE[idx_ge1] *
                                                  log2(Results_combined$RE[idx_ge1])) / Results_combined$RE[idx_ge1]
  
  ## ---- reorder columns: put lower/upper SEs before letters ----
  # find current column positions
  # We'll place Lower.se.RE, Upper.se.RE, Lower.se.log2FC, Upper.se.log2FC before 'letters'
  cols <- colnames(Results_combined)
  letters_pos <- match("letters", cols)
  new_order <- c(cols[1:(letters_pos - 1)],
                 "Lower.se.RE", "Upper.se.RE", "Lower.se.log2FC", "Upper.se.log2FC",
                 cols[letters_pos:length(cols)])
  # deduplicate and keep only existing columns
  new_order <- unique(new_order[new_order %in% cols])
  Results_final <- Results_combined[, new_order, drop = FALSE]
  Results_final$RowNames <- NULL
  
  ## ---- prepare output ----
  # remove the temporary column T from original data for output (you created x$T earlier)
  xx <- x[, setdiff(names(x), "T"), drop = FALSE]
  
  outlist2 <- structure(list(Final_data = xx,
                             lmCRD = lm_fit,
                             ANOVA = anovaCRD,
                             Results = Results_final),
                        class = "XX")
  
  print.XX <- function(outlist2) {
    print(outlist2$ANOVA)
    cat("\n", sep = '', "Relative expression (DCt method)", "\n")
    print(outlist2$Results)
    invisible(outlist2)
  }
  
  print.XX(outlist2)
}
