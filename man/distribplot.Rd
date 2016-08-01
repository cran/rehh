\name{distribplot}
\alias{distribplot}
\title{Distribution of standardized iHS or Rsb values}
\description{Plot the observed distribution of standardized iHS or Rsb values together 
             with the expected standard Gaussian distribution}
\usage{
distribplot(data,lty = 1,lwd = 1.5,col = c("blue","red"),
main = "Genome-wide distribution",xlab = "",cex.main = 1.5,cex.lab = 1.25,qqplot = TRUE)
}
\arguments{
  \item{data}{A vector of iHS, Rsb or XPEHH values.}
  \item{col}{A vector describing color of the Observed and expected Gaussian distribution}
  \item{main}{Character string for the plot legend}
  \item{xlab}{Character string for the X-axis legend}
  \item{cex.lab}{Size of axis legends}
  \item{cex.main}{Size of the main legend}
  \item{lty}{Line Type}
  \item{lwd}{Line Width}
  \item{qqplot}{If TRUE a qq-plot is drawn}
}
\value{The function returns a plot.}
\seealso{scan_hh,ihh2ihs,ies2rsb,ihsplot,rsbplot,ies2xpehh,xpehhplot}
\examples{
data(wgscan.cgu) 
## results from a genome scan (44,057 SNPs) see ?wgscan.eut for details
val.ihs<-ihh2ihs(wgscan.cgu)$iHS[,3]
##standardize
distribplot(val.ihs,main="iHS (CGU population)")
dev.off()
}
