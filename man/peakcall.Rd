% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peakcall.R
\name{peakcall}
\alias{peakcall}
\title{Peak calling}
\usage{
peakcall(res, flank = 1, fdrcut = 1e-05)
}
\arguments{
\item{res}{Result matrix returned by SCATE}

\item{flank}{Numeric variable of the number of flanking bins for each bin. For each bin, an averaged signal of itself and the flanking bins will be calculated and compared to a background distribution.}

\item{fdrcut}{Numeric variable of FDR cutoff. Bins passing the FDR cutoff will be peaks.}
}
\value{
A list with length equal to the number of clusters. Each element is a data frame with five columns: chromosome name, starting location, ending location, FDR and signal. The data frame is ordered by FDR then by signal.
}
\description{
Peak calling function
}
\details{
This function performs peak calling for signal generated by SCATE
}
\examples{
peakcall(SCATE(GRanges(seqnames="chr1",IRanges(start=1:100+1e6,end=1:100+1e6)),clunum=5000,genome='mm10'))
}
\author{
Zhicheng Ji, Weiqiang Zhou, Hongkai Ji <zji4@zji4.edu>
}
