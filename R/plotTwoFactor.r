#' @title Bar plot of gene expression from an expression table of a two-factorial experiment data
#' 
#' @description Bar plot of the relative expression (\eqn{\Delta C_T} method) of a gene along with the standard error (se), 95\% confidence interval (ci) and significance
#' 
#' @details The \code{plotTwoFactor} function generates the bar plot of the average fold change for target genes along with the significance, standard error (se) and the 95\% confidence interval (ci) as error bars.
#' 
#' @author Ghader Mirzaghaderi
#' @export plotTwoFactor
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @param data A data.frame such as the expression result of \code{ANOVA_DDCt(x)} or \code{ANOVA_DCt(x)}, etc. functions.
#' @param x_col column number of the x-axis factor.
#' @param y_col column number of the y-axis factor.
#' @param group_col column number of the grouping factor.
#' @param Lower.se_col The column number of the data.frame used for the lower error bar.
#' @param Upper.se_col The column number of the data.frame used for the upper error bar.
#' @param letters_col The column number of the data.frame used as the result of statistical comparing and grouping.
#' @param show.groupingLetters a logical variable. If TRUE, mean grouping letters (the results of statistical comparison) are added to the bars.
#' @return Bar plot of the average fold change for target genes along with the standard error or 95\% confidence interval as error bars.
#' @examples
#' 
#' a <- ANOVA_DCt(data_2factorBlock, block = "Block", numberOfrefGenes = 1)
#' data <- a$Results
#' 
#' plotTwoFactor(
#'   data = data,
#'   x_col = 2,
#'   y_col = 3,
#'   group_col = 1,
#'   Lower.se_col = 8,
#'   Upper.se_col = 9,
#'   letters_col = 12,
#'   show.groupingLetters = TRUE)
#'               
#'               
#'               
#' # Combining FC results of two different genes:
#' a <- REPEATED_DDCt(data_repeated_measure_1,
#'                    numberOfrefGenes = 1,
#'                    factor = "time", block = NULL, plot = FALSE)
#' 
#' b <- REPEATED_DDCt(data_repeated_measure_2,
#'                    factor = "time",
#'                    numberOfrefGenes = 1, block = NULL, plot = FALSE)
#' 
#' a1 <- a$FC_statistics_of_the_main_factor
#' b1 <- b$FC_statistics_of_the_main_factor
#' c <- rbind(a1, b1)
#' c$gene <- factor(c(1,1,1,2,2,2))
#' c
#' 
#' plotTwoFactor(
#'   data = c,
#'   x_col = 1,
#'   y_col = 2,
#'   group_col = 13,
#'   Lower.se_col = 9,
#'   Upper.se_col = 10,
#'   letters_col = 5,
#'   show.groupingLetters = TRUE)
#'
#'


plotTwoFactor <- function(data, 
                          x_col,        # x-axis factor
                          y_col,        # bar height
                          group_col,    # fill groups
                          Lower.se_col, # lower SE column
                          Upper.se_col, # upper SE column
                          letters_col = NULL,
                          show.groupingLetters = TRUE){
  
  # Extract Column Names by Index
  x_name      <- names(data)[x_col]
  y_name      <- names(data)[y_col]
  group_name  <- names(data)[group_col]
  lower  <- names(data)[Lower.se_col]
  upper  <- names(data)[Upper.se_col]
  letter_name <- if (!is.null(letters_col)) names(data)[letters_col] else NULL
  
  # compute ymin and ymax BEFORE ggplot
  data$ymin <- data[[lower]]
  data$ymax <- data[[upper]]

  
  #Base ggplot
  p <- ggplot(data,
              aes(x = .data[[x_name]],
                  y = .data[[y_name]],
                  fill = .data[[group_name]])) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.8)
  

  
  
  #Error Bars
  p <- p + geom_errorbar(
    aes(ymin = ymin,
        ymax = ymax),
    width = 0.15,
    position = position_dodge(width = 0.8))
  
  #Add Grouping Letters if Provided
  if (show.groupingLetters & !is.null(letters_col)) {

    letters_name <- names(data)[letters_col]
    p <- p +
      geom_text(
        aes(
          label = .data[[letters_name]],
          y = ifelse(
            .data[[y_name]] < 0,
            ymin - 0.3,   # negative bars
            ymax + 0.3   # positive bar
          )
        ), position = position_dodge(width = 0.8)
      )
    
    
  }
  
  return(p)
}

