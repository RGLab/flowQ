\name{qaProcess.timeline}
\alias{qaProcess.timeline}

\title{  Create QA process of type 'timeline' }

\description{
  
  This function takes a \code{\link[flowCore:flowSet-class]{flowSet}} as
  input and 
  creates all necessary output for a 'timeline' type QA process. Objects
  created by this function can be laid out as HTML using
  \code{\link[flowQ]{writeQAReport}}.
  
}

\usage{

qaProcess.timeline(set, channels=NULL, outdir, cutoff=1,
  name="time line", sum.dimensions=NULL, det.dimensions=c(3,2),
  pdf=TRUE, \dots)

}


\arguments{
  
  \item{set}{ A \code{\link[flowCore:flowSet-class]{flowSet}} }
  
  \item{channels}{ A character vector of channel names for which the
    qaReport is to be produced. Will use all available channels except
    for the time channel if \code{NULL}.}
  
  \item{outdir}{ The directory to which the graphical output is to be
    saved. If multiple QA processes are to be combined, make sure to use
    the same directory for all of them.}

  \item{cutoff}{ The threshold at which the QA criterion is considered
    to be failed. An absolute value in the timeline deviation score as
    computed by the
    \code{\link[flowViz:timeLinePlot-methods]{timeLinePlot}} function.} 
  
  \item{name}{ The name of the process used for the headings in the HTML
    output.}
  
  \item{sum.dimensions, det.dimensions}{ The dimensions of the pdf device
    in inches used for the summary and the detailed plots.}
  
  \item{pdf}{ Logical indicating whether to create vectorized versions of
     images for this quality process. This should be set to \code{FALSE}
     if disk space is critical, since the pdf versions of the image
     consume much more space on the hard drive compared to the bitmap
     version.}
  
  \item{\dots}{ Further arguments.}

 }
 \details{
   
  QA processes of type 'timeline' detect unusal patterns in the
  acquisition of fluorescense and light scatter measurements over
  time. These are detected dynamically by identifying trends in the
  signal intensity over time or local changes in the measurement
  intensities. The underlying hypothesis is that measurement values are
  acquired randomly, hence there shouldn't be any correlation to time.
  
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
  \code{\link[flowQ]{qaProcess.timeflow}},
  \code{\link[flowQ]{qaProcess.cellnumber}}

}


\examples{
\dontrun{
data(GvHD)
GvHD <- transform(GvHD, "FL1-H"=asinh(`FL1-H`), "FL2-H"=asinh(`FL2-H`))
dest <- tempdir()
qp <- qaProcess.timeline(GvHD, channel="FL1-H", outdir=dest, cutoff=1)
qp
}

}

\keyword{ dynamic }


