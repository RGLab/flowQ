\name{qaProcess.cellnumber}
\alias{qaProcess.cellnumber}

\title{  Create QA process of type 'cellnumber' }

\description{
  
  This function takes a \code{\link[flowCore:flowSet-class]{flowSet}} as
  input and 
  creates all necessary output for a 'cellnumber' type QA
  process. Objects created by this function can be laid out as HTML
  using \code{\link[flowQ]{writeQAReport}}.
  
}

\usage{

qaProcess.cellnumber(set, grouping=NULL, outdir, cFactor=2,
                     absolute.value=NULL, two.sided = FALSE,
                     name="cell number", sum.dimensions=NULL,
                     pdf=TRUE, \dots)

}

\arguments{
  
  \item{set}{ A \code{\link[flowCore:flowSet-class]{flowSet}}. }

  \item{grouping}{ A character vector defining one of the variables in
  the phenoData of \code{set} used as a grouping variable. If this
  argument is used and if \code{absolute.value} is \code{NULL}, outlier
  detection will be performed within groups rather than across all
  samples.}
  
  \item{outdir}{ The directory to which the graphical output is to be
    saved. If multiple QA processes are to be combined, make sure to use
    the same directory for all of them.}

  \item{cFactor}{ The outlier threshold at which the QA criterion is
    considered to have failed. This is essentially the factor of
    standard deviations away from the average number of cells per
    sample, either in both directions if \code{two.sided=TRUE} or only
    towards smaller event numbers if \code{two.sided=FALSE}.}

  \item{absolute.value}{ An absolute event count below which the QA
    criterion is considered to be failed. If this argument is not
    \code{NULL}, \code{cFactor} and \code{two.sided} are ignored.}

  \item{two.sided}{ Perform a two-sided outlier detection, i.e., report
    both unproportionally high and low event numbers.}
  
  \item{name}{ The name of the process used for the headings in the HTML
    output.}
  
   \item{sum.dimensions}{ The dimensions of the pdf deviced in inches
     used for the summary plot. }
   
   \item{pdf}{ Logical indicating whether to create vectorized versions
     of images for this quality process. This should be set to
     \code{FALSE} if disk space is critical, since the pdf versions of
     the image consume much more space on the hard drive compared to the
     bitmap version.}
   
   \item{\dots}{Further arguments.}
   
 }
 \details{
   
  QA processes of type 'cellnumber' detect aberations in the number of
  events per sample. These are either determined dynamically as outliers
  from the typical distribution of event counts for the whole set, or,
  if \code{absolute.value} is not \code{NULL}, by an absolute cutoff
  value. If there is a natural grouping among the samples, this can be
  specified using the \code{grouping} argument. In this case, the
  outlier detection will be performed within its respective group for a
  particular sample.
  
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
  \code{\link[flowQ]{qaProcess.marginevents}},
  \code{\link[flowQ]{qaProcess.timeflow}}, 
  \code{\link[flowQ]{qaProcess.timeline}}

}


\examples{
\dontrun{
data(GvHD)
dest <- file.path(tempdir(), "flowQ")
qp <- qaProcess.cellnumber(GvHD, outdir=dest, cFactor=2)
qp
}

}

\keyword{ dynamic }


