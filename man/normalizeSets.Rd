\name{normailzeSets}
\alias{normalizeSets}

\title{Locate Duplicated Parameters }
\description{
Aligns up peaks for parameters across multiple flowSets using the
\code{\link[flowStats]{warpSet}} method from the \code{flowStats} package.

}
\usage{
  normalizeSets(flowList,dupes,peaks=NULL)
}

\arguments{
  \item{flowList}{ A list of \code{flowSet} }
  \item{dupes}{ Stains that are to be aligned across multiple \code{flowSet} }
  \item{peaks}{ Number of peaks expected in the distribution. 
   If NULL is passed as input, The number of peaks to align are
  estimated from the distribution } 
}
 
\details{  
  This function can be used to normalize dye information across multiple
  flowSets so that the peaks in the distributions align.  The method
  makes use of information provided in the description field of each
  flowFrame to align up floursecense  nformation from multiple \code{flowSet}.
}

\value{
  
  A \code{list} of objects of class
  \code{\link[flowCore:flowSet-class]{flowSet}}. 
}

\author{ Nishant Gopalakrishnan }


\seealso{
  \code{\link[flowQ]{locateDuplicatedParameters}}
  \code{\link[flowStats]{warpSet}}
}


\examples{

\dontrun{


}

}

\keyword{ dynamic }


 
