% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qpcrTTESTplot.r
\name{qpcrTTESTplot}
\alias{qpcrTTESTplot}
\title{Bar plot of the average fold change (\eqn{\Delta \Delta C_T} method) of target genes}
\usage{
qpcrTTESTplot(
  x,
  order = "none",
  numberOfrefGenes,
  paired = FALSE,
  var.equal = TRUE,
  p.adj = "BH",
  width = 0.5,
  fill = "skyblue",
  y.axis.adjust = 0,
  y.axis.by = 2,
  letter.position.adjust = 0.3,
  ylab = "Average Fold Change",
  xlab = "none",
  fontsize = 12,
  fontsizePvalue = 7,
  axis.text.x.angle = 0,
  axis.text.x.hjust = 0.5,
  errorbar = "se"
)
}
\arguments{
\item{x}{a data frame. The data frame consists of 4 columns belonging to condition levels, 
E (efficiency), genes and Ct values, respectively. Each Ct in the following data frame is 
the mean of technical replicates. Complete amplification efficiencies of 2 is assumed here 
for all wells but the calculated efficienies can be used we well. We use this data set for 
fold change expression analysis of the target genes in treatment condition compared to 
normal condition. See \href{../doc/vignette.html}{\code{vignette}} for details.}

\item{order}{a vector determining genes order on the output graph.}

\item{numberOfrefGenes}{number of reference genes. Up to two reference genes can be handled.}

\item{paired}{a logical indicating whether you want a paired t-test.}

\item{var.equal}{a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.}

\item{p.adj}{Method ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none") for adjusting p values.}

\item{width}{a positive number determining bar width.}

\item{fill}{specify the fill color for the columns of the bar plot.}

\item{y.axis.adjust}{a negative or positive value for reducing or increasing the length of the y axis.}

\item{y.axis.by}{determines y axis step length}

\item{letter.position.adjust}{adjust the distance between the signs and the error bars.}

\item{ylab}{the title of the y axis}

\item{xlab}{the title of the x axis}

\item{fontsize}{fonts size of the plot}

\item{fontsizePvalue}{font size of the pvalue labels}

\item{axis.text.x.angle}{angle of x axis text}

\item{axis.text.x.hjust}{horizontal justification of x axis text}

\item{errorbar}{Type of error bar, can be \code{se} or \code{ci}.}
}
\value{
Bar  plot of the average fold change for target genes along with the significance and the 95 percent CI as error bars.
}
\description{
Bar plot of the fold change (\eqn{\Delta \Delta C_T} method) values for for any number of target 
genes under a two-level conditional experimental (e.g. control and treatment).
}
\details{
The \code{qpcrTTESTplot} function applies a t.test based analysis to any number of target genes 
along with one or two reference gene(s), that have been evaluated under control and treatment conditions. 
It returns the bar plot of the fold change (FC) values for target genes along with the 95\% CI and significance. 
Sampling may be unpaired or paired. Unpaired and paired samples are commonly analyzed using unpaired and 
paired t-test, respectively.Paired samples in quantitative PCR refer to two sample data that are collected 
from one set of individuals 
at two different conditions, for example before and after a treatment or at two different time points. While 
for unpaired samples, two sets of individuals are used: one under untreated and the other set under treated 
condition. Paired samples allow to compare gene expression changes within the same individual, reducing 
inter-individual variability.
}
\examples{

# See a sample data frame
data_ttest


qpcrTTESTplot(data_ttest, 
              numberOfrefGenes = 1,
              errorbar = "ci",
              p.adj = "BH")


# Producing the plot
qpcrTTESTplot(data_ttest,
              numberOfrefGenes = 1,
              order = c("C2H2-01", "C2H2-12", "C2H2-26"),
              paired = FALSE,
              var.equal = TRUE,
              width = 0.5,
              fill = "skyblue",
              y.axis.adjust = 0,
              y.axis.by = 2,
              letter.position.adjust = 0.3,
              ylab = "Fold Change in Treatment vs Control",
              xlab = "Gene",
              errorbar = "se")


}
\author{
Ghader Mirzaghaderi
}
