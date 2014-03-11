\name{evaluateProcess}
\alias{evaluateProcess}

\title{ Evaluate QA processes }
\description{
  Re-evaluate an object of class \code{\link{qaProcess}}, e.g. for the
  case that a threshold value has changed.
}
\usage{
evaluateProcess(process, thresh, ...)
}
%
\arguments{
  \item{process}{ An object of class \code{\link{qaProcess}}. }
  \item{thresh}{ The new treshold on which the process is to be evaluated. }
  \item{\dots}{ Further arguments that are passed on to the individual
    functions for each QA process type.}
}
\details{
  It is sometimes useful to update the state of aggregators in a
  \code{\link{qaProcess}} , for instances after changing the threshold
  value, without having to recompute all images, which
  can be very time consuming.
}
\value{
 An updated object of class \code{\link[flowQ:qaProcess-class]{qaProcess}}
}

\author{ Florian Hahne }
\note{
  This function needs to be extended for new types of
  \code{\link[flowQ:qaProcess-class]{qaProcess}}. 
}
\seealso{ \code{\link[flowQ:qaProcess-class]{qaProcess}},
  \code{\link[flowQ]{writeQAReport}}}
\examples{
\dontrun{
data(GvHD)
dest <- tempdir()
qp1 <- qaProcess.timeline(GvHD[1:3], channel="FL1-H", outdir=dest,
         cutoff=1)
evaluateProcess(qp1, thresh=4)
}
}

\keyword{ dynamic }
