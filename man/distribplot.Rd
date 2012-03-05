\name{distribplot}
\alias{distribplot}
\title{Distribution of standardized iHS or Rsb values}
\description{Plot the observed distribution of standardized iHS or Rsb values together with the expected standard Gaussian distribution}
\usage{
distribplot(data,col=c("blue","red"),
main="iHS distribution",xlab="iHS")
}
\arguments{
  \item{data}{A vector of iHS or Rsb values.}
  \item{col}{A vector describing color of the Observed and expected Gaussian distribution}
  \item{main}{Character string for the plot legend}
  \item{xlab}{Character string for the X-axis legend}
}
\value{The function returns a plot.}
\seealso{scan_hh,ihh2ihs,ies2rsb,ihsplot,rsbplot}
\examples{
data(wgscan.cgu) 
## results from a genome scan (44,057 SNPs) see ?wgscan.eut for details
val.ihs<-ihh2ihs(wgscan.cgu)$res.ihs[,3]
##standardize
distribplot(val.ihs,main="iHS (CGU population)")
dev.off()
}
