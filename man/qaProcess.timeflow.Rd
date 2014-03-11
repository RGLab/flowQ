\name{qaProcess.timeflow}
\alias{qaProcess.timeflow}

\title{  Create QA process of type 'timeflow' }

\description{

  This function takes a \code{\link[flowCore:flowSet-class]{flowSet}} as
  input and 
  creates all necessary output for a 'timeflow' type QA
  process. Objects created by this function can be laid out as HTML
  using \code{\link[flowQ]{writeQAReport}}.
  
}

\usage{

qaProcess.timeflow(set, outdir, cutoff=2, name="time flow",
sum.dimensions=NULL, det.dimensions=c(3,2), pdf=TRUE, \dots)

}

\arguments{
  
  \item{set}{ A \code{\link[flowCore:flowSet-class]{flowSet}}.}
  
  \item{outdir}{ The directory to which the graphical output is to be
    saved. If multiple QA processes are to be combined, make sure to use
    the same directory for all of them.}

  \item{cutoff}{ The threshold at which the QA criterion is considered
    to have failed. An absolute value in the timeline deviation score as
    computed by the
    \code{\link[flowViz:timeLinePlot-methods]{timeLinePlot}} function.} 
  
  \item{name}{ The name of the process used for the headings in the HTML
    output.}
  
   \item{sum.dimensions, det.dimensions}{ The dimensions of the pdf
     device in inches used for the summary and the detailed plots.}
   
   \item{pdf}{ Logical indicating whether to create vectorized versions
     of images for this quality process. This should be set to
     \code{FALSE} if disk space is critical, since the pdf versions of
     the image consume much more space on the hard drive compared to the
     bitmap version.}
   
   \item{\dots}{Further arguments.}

 }
 \details{
   
  QA processes of type 'timeflow' detect disturbances in the flow of
  cells over time. These indicate problems with the cell suspension,
  clogging of the needle or similar issues. If the flow rate is to high,
  the frequency of measuring cell douplettes increases.
  
  For more details on how to layout
  \code{\link[flowQ:qaProcess-class]{qaProcess}} objects to 
  HTML, see \code{\link[flowQ]{writeQAReport}} and
  \code{\link[flowQ]{qaReport}}. 
  
}
\value{
  
  An object of class \code{\link[flowQ:qaProcess-class]{qaProcess}}.
  
}

\author{ Florian Hahne }



\seealso{

  \code{\link[flowQ]{writeQAReport}},
  \code{\link[flowQ]{qaReport}},
  \code{\link[flowQ:qaProcess-class]{qaProcess}},
  \code{\link[flowQ]{qaProcess.marginevents}}
  \code{\link[flowQ]{qaProcess.timeline}},
  \code{\link[flowQ]{qaProcess.cellnumber}}
  
}


\examples{

\dontrun{
data(GvHD)
dest <- tempdir()
qp <- qaProcess.timeflow(GvHD, outdir=dest, cutoff=1)
qp
}

}

\keyword{ dynamic }


