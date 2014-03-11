\name{locateDuplicatedParameters}
\alias{locateDuplicatedParameters}

\title{Locate Duplicated Parameters }
\description{
Identifies dyes that are repeated arcross multiple flowSets
using information provided in the description field of each flowFrame.
}

\usage{
  locateDuplicatedParameters(flowList)
}

\arguments{
  \item{flowList}{ A list of \code{\link[flowCore:flowSet-class]{flowSet}} }
}
 \details{
  This function is used internally by the qaProcess functions to identify dyes
  that are repeated across multiple aliquots. The function takes in  a
  list of \code{flowSet} and identifies stains that are repeated across
  multiple flowSets (aliquots). These parameters are expected to have
  the same distribution across aliquots and can be used for Quality
  control. The method looks for the information provided in the
  description field of each flowFrame. 
  
}
\value{
  
  An object of class \code{\link[flowQ:qaProcess-class]{qaProcess}}.  
}

\author{ Nishant Gopalakrishnan }


\seealso{
  \code{\link[flowQ]{normalizeSets}}
}


\keyword{ dynamic }


 
