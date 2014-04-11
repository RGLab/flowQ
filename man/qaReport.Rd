\name{qaReport}
\alias{qaReport}

\title{ Create HTML report using one or several QA process function(s) }
\description{
  
  This function combines all graphical output of multiple QA process
  functions for one \code{\link[flowCore:flowSet-class]{flowSet}} in a single
  hyperlinked HTML document.

}
\usage{
qaReport(set, qaFunctions, outdir = "./qaReport",
         argLists, grouping = NULL, \dots)
}

\arguments{
  \item{set}{ A \code{\link[flowCore:flowSet-class]{flowSet}} }
  
  \item{qaFunctions}{ A character vector of the names of QA process
    functions to be used }
  
  \item{outdir}{ The directory to which the HTML report is to be
    saved. }
  
  \item{argLists}{ lists of argument lists for each of the QA process
    functions specified via \code{qaFunctions} }
  
  \item{grouping}{ A character scalar indicating a variable in the
    \code{\link[flowCore:flowSet-class]{flowSet}}'s \code{phenoData}
    that is used as a grouping factor in the output. }

  \item{\dots}{Further arguments that are passed on to
    \code{\link[flowQ]{writeQAReport}}.}
  
}
\details{
  
  This is a simple convenience function to produce HTML QA reports
  for a single \code{\link[flowCore:flowSet-class]{flowSet}} given a list of QA
  process functions. For more fine-grained control use function 
  \code{\link[flowQ]{writeQAReport}} directly.

  An entry point to the output of this function can be found at
  \code{outdir/index.html}.
}
\value{
 The function is called for it's side effects
}

\author{ Florian Hahne }

\seealso{ \code{\link[flowQ]{qaProcess.marginevents}},
  \code{\link[flowQ]{writeQAReport}},
  \code{\link[flowQ:qaProcess-class]{qaProcess}},
  \code{\link[flowQ]{qaProcess.timeline}}}
\examples{
\dontrun{
data(GvHD)
GvHD <- transform(GvHD, "FL1-H"=asinh(`FL1-H`), "FL2-H"=asinh(`FL2-H`))
dest <- tempdir()
qaReport(GvHD, c("qaProcess.timeline", "qaProcess.marginevents"), dest,
 list(list(channel="FL1-H", cutoff=1), list(channels=c("FL1-H",
  "FL2-H"), cFactor=4)))
browseURL(file.path(dest, "index.html"))
}
}

\keyword{ IO }
\keyword{ dynamic }
