#' @title Multiple plot function
#' @description \code{multiplot} function combines multiple ggplot objects into a single plate.
#' @details Combining multiple ggplot objects into a single plate.
#' @author gist.github.com/pedroj/ffe89c67282f82c1813d
#' @export multiplot
#' @import grid
#' @param ... ggplot objects can be passed in ... or to plotlist (as a list of ggplot objects)
#' @param cols Number of columns in the panel
#' @return A  multiple-plots plate
#' @examples
#' 
#' a <- TTEST_DDCt(data_ttest, 
#'         numberOfrefGenes = 1)
#' p1 <- a$plot
#' 
#' out2 <- ANOVA_DCt(data_1factor, numberOfrefGenes = 1, block = NULL)$Result
#' p2 <- plotOneFactor(out2,
#'         x_col= 1,
#'         y_col= 2,
#'         Lower.se_col = 7,
#'         Upper.se_col = 8,
#'         letters_col = 11,
#'         show.groupingLetters = TRUE)
#'                     
#' multiplot(p1, p2, cols=2)
#' 
#' multiplot(p1, p2, cols=1)
#'
#'



multiplot <- function(..., cols=1) {
  
  # Make a list from the ... arguments
  plots <- c(list(...))
  
  numPlots = length(plots)
  
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow = TRUE)
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
