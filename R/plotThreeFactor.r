#' @title Bar plot of gene expression from an expression table of a three-factorial experiment data
#' 
#' @description Bar plot of the relative expression (\eqn{\Delta C_T} method) of a gene along with the confidence interval and significance
#' 
#' @details The \code{plotThreeFactor} function generates the bar plot of the average fold change for target genes along with the significance, standard error (se) and the 95\% confidence interval (ci).
#' 
#' @author Ghader Mirzaghaderi
#' @export plotThreeFactor
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @param data A data.frame such as the expression result of \code{ANOVA_DDCt(x)} or \code{ANOVA_DCt(x)}, etc. functions.
#' @param x_col column number of the x-axis factor.
#' @param y_col column number of the y-axis factor.
#' @param group_col column number of the grouping factor.
#' @param facet_col column number of the faceting factor.
#' @param Lower.se_col The column number of the data.frame used for the lower error bar.
#' @param Upper.se_col The column number of the data.frame used for the upper error bar.
#' @param letters_col The column number of the data.frame used as the result of statistical comparing and grouping.
#' @param show.groupingLetters a logical variable. If TRUE, mean grouping letters (the results of statistical comparison) are added to the bars.
#' @return Bar plot of the average fold change for target genes along with the standard error or 95\% confidence interval as error bars.
#' @examples
#' 
#'
#'
#' res <- ANOVA_DCt(data_3factor, 
#'       numberOfrefGenes = 1, 
#'       block = NULL)
#'       
#' data <- res$Results
#' 
#' plotThreeFactor(data, 
#'                 3,   # x-axis factor
#'                 5,   # bar height
#'                 1,   # fill groups
#'                 2,   # facet grid
#'                 11,   # lower SE column
#'                 12,   # upper SE column
#'                 letters_col = 13,
#'                 show.groupingLetters = TRUE)
#'
#'
#'


plotThreeFactor <- function(data, 
                            x_col,        # x-axis factor
                            y_col,        # bar height
                            group_col,    # fill groups
                            facet_col,    # facet grid
                            Lower.se_col, # lower SE column
                            Upper.se_col, # upper SE column
                            letters_col = NULL,
                            show.groupingLetters = TRUE) {
  
  # Extract column names
  x_name     <- names(data)[x_col]
  y_name     <- names(data)[y_col]
  group_name <- names(data)[group_col]
  facet_name <- names(data)[facet_col]
  lower_name <- names(data)[Lower.se_col]
  upper_name <- names(data)[Upper.se_col]
  
  # compute ymin and ymax BEFORE ggplot
  data$ymin <- data[[lower_name]]
  data$ymax <- data[[upper_name]]


  
  # Base plot
  p <- ggplot(data,
              aes(x = .data[[x_name]],
                  y = .data[[y_name]],
                  fill = .data[[group_name]])) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.15,
                  position = position_dodge(width = 0.9)) +
    facet_wrap(vars(.data[[facet_name]])) +
    theme_bw() +
    labs(x = x_name,
         y = y_name,
         fill = group_name)
  
  # ---- Optional grouping letters ----
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
            ymin - 0.5,   # negative bars
            ymax + 0.3   # positive bar
          )
        ),
        position = position_dodge(width = 0.9),
        vjust = 0
      )
  }
  
  return(p)
}
