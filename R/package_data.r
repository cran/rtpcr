#' Sample data (amplification efficiency)
#'
#' A sample qPCR dataset for demonstrating efficiency calculation.
#'
#' @format A data frame with 21 observations and 3 variables:
#' \describe{
#'   \item{dilutions}{Dilution factor}
#'   \item{C2H2.26}{Target gene}
#'   \item{GAPDH}{Reference gene}
#' }
#'
#' @source Where the data comes from (if applicable)
"data_efficiency"

#' Sample data (one factor three levels)
#'
#' A sample dataset for demonstration purposes.
#'
#' @format A data frame with 9 observations and 6 variables:
#' \describe{
#'   \item{SA}{First experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
"data_1factor"

#' Sample data (two factor)
#'
#' A sample dataset for demonstration purposes.
#'
#' @format A data frame with 18 observations and 7 variables:
#' \describe{
#'   \item{Genotype}{First experimental factor}
#'   \item{Drought}{Second experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
"data_2factor"

#' Sample data (two factor with blocking factor)
#'
#' A sample qPCR data set with blocking factor.
#'
#' @format A data frame with 18 observations and 8 variables:
#' \describe{
#'   \item{factor1}{First experimental factor}
#'   \item{factor2}{Second experimental factor}
#'   \item{block}{Second experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
"data_2factorBlock"

#' Sample data (three factor)
#'
#' A sample dataset for demonstration purposes.
#'
#' @format A data frame with 36 observations and 8 variables:
#' \describe{
#'   \item{Genotype}{First experimental factor}
#'   \item{Drought}{Second experimental factor}
#'   \item{SA}{Third experimental factor}
#'   \item{Rep}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source Not applicable
"data_3factor_a"

#' Sample data (three factor)
#'
#' A sample real time PCR data for demonstration purposes.
#'
#' @format A data frame with 36 observations and 8 variables:
#' \describe{
#'   \item{Type}{The first experimental factor}
#'   \item{Conc}{The second experimental factor}
#'   \item{SA}{The third experimental factor}
#'   \item{Replicate}{Biological replicates}
#'   \item{EPO}{Mean amplification efficiency of PO gene}
#'   \item{POCt}{Ct values of PO gene. Each is the mean of technical replicates}
#'   \item{EGAPDH}{Mean amplification efficiency of GAPDH gene}
#'   \item{GAPDHCt}{Ct values of GAPDH gene. Each is the mean of technical replicates}
#' }
#'
#' @source University of Kurdistan
"data_3factor_b"

#' Sample data (one factor-two level qPCR)
#'
#' A sample data for demonstrating qPCR data analysis.
#'
#' @format A data frame with 24 observations and 4 variables:
#' \describe{
#'   \item{Condition}{Experimental conditions}
#'   \item{Gene}{Genes}
#'   \item{E}{Amplification efficiency}
#'   \item{Ct}{Ct values}
#' }
#'
#' @source University of Kurdistan
"data_ttest"

#' Sample data (one target, two reference)
#'
#' One target and two reference gens for demonstrating qPCR data analysis.
#'
#' @format A data frame with 18 observations and 4 variables:
#' \describe{
#'   \item{Condition}{Experimental conditions}
#'   \item{Gene}{Genes}
#'   \item{E}{Amplification efficiency}
#'   \item{Ct}{Ct values}
#' }
#'
#' @source Not applicable
"data_ttest2"

#' Sample data (with technical replicates)
#'
#' A sample data for calculating biological replicated.
#'
#' @format A data frame with 18 observations and 9 variables:
#' \describe{
#'   \item{factor1}{experimental factor}
#'   \item{factor2}{experimental factor}
#'   \item{factor3}{experimental factor}
#'   \item{biolrep}{biological replicate}
#'   \item{techrep}{technical replicates}
#'   \item{Etarget}{Amplification efficiency of target gene}
#'   \item{targetCt}{Ct of target gene}
#'   \item{Eref}{Amplification efficiency of reference gene}
#'   \item{refCt}{Ct of reference gene}
#' } 
#'
#' @source Not applicable
"data_withTechRep"


#' Sample data (with technical replicates)
#'
#' A sample data for calculating biological replicated.
#'
#' @format A data frame with 72 observations and 8 variables:
#' \describe{
#'   \item{factor1}{experimental factor}
#'   \item{DS}{DS}
#'   \item{biolRep}{biological replicate}
#'   \item{techRep}{technical replicates}
#'   \item{APOE_efficiency}{Amplification efficiency of APOE gene}
#'   \item{APOE_Ct}{Ct of APOE gene}
#'   \item{GAPDH_efficiency}{Amplification efficiency of GAPDH gene}
#'   \item{GAPDH_Ct}{Ct of GAPDH gene}
#' }
#'
#' @source Lee et al, (2020) <doi:10.12688/f1000research.23580.2>
"Lee_etal2020qPCR"
