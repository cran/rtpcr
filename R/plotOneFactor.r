#' @title Bar plot of gene expression from an expression table of a single-factor experiment data
#' 
#' @description Bar plot of the relative expression of a gene along with the standard error (se), 95\% confidence interval (ci) and significance. 
#' 
#' @details The \code{plotOneFactor} function generates the bar plot of the  fold change for target genes along with the significance and the 95\% confidence interval as error bars.
#' 
#' @author Ghader Mirzaghaderi
#' @export plotOneFactor
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @param data A data.frame such as the expression result of \code{ANOVA_DDCt(x)} or \code{ANOVA_DCt(x)}, etc. functions.
#' @param show.groupingLetters a logical variable. If TRUE, mean grouping letters (the results of statistical comparison) are added to the bars.
#' @param x_col The column number of the data.frame used for x axis.
#' @param y_col The column number of the data.frame used for y axis.
#' @param Lower.se_col The column number of the data.frame used for the lower error bar.
#' @param Upper.se_col The column number of the data.frame used for the upper error bar.
#' @param letters_col The column number of the data.frame used as the result of statistical comparing and grouping.
#' @return Bar plot of the average fold change for target genes along with the significance and the standard error or 95\% confidence interval as error bars.
#' @examples
#'
#' # Before plotting, the result needs to be extracted as below:
#' res <- ANOVA_DCt(data_1factor, numberOfrefGenes = 1, block = NULL)$Result
#'
#' # Bar plot
#' plotOneFactor(res, 1, 2, 7, 8, 11,
#'     show.groupingLetters = TRUE)
#'
#'


plotOneFactor <- function(data, 
                          x_col, 
                          y_col,
                          Lower.se_col,
                          Upper.se_col,
                          letters_col = NULL,
                          show.groupingLetters = TRUE) {

  
# column names from index
  x_name  <- names(data)[x_col]
  y_name  <- names(data)[y_col]
  lower   <- names(data)[Lower.se_col]
  upper   <- names(data)[Upper.se_col]
  
# compute ymin and ymax BEFORE ggplot
  data$ymin <- ifelse(
    data[[y_name]] < 0,
    data[[lower]],   # negative bars
    data[[lower]]    # positive bars
  )
  
  data$ymax <- ifelse(
    data[[y_name]] < 0,
    data[[upper]],   # negative bars
    data[[upper]]    # positive bars
  )
  
  # Base plot
  p <- ggplot(data, aes(x = .data[[x_name]], y = .data[[y_name]])) +
    geom_col() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1)
  
  # Add grouping letters
  if (show.groupingLetters) {
    if (is.null(letters_col)) {
      stop("letters_col must be provided when show.groupingLetters = TRUE")
    }
    
    letters_name <- names(data)[letters_col]
    
    p <- p +
      geom_text(
        aes(
          label = .data[[letters_name]],
          y = ifelse(
            .data[[y_name]] < 0,
            .data[[lower]] - 0.2,   # negative bars
            .data[[upper]] + 0.2   # positive bar
          )
        )
      )
  }
  
  return(p)
}
