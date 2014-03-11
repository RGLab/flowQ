\name{preProcessFlowList}
\alias{preProcessFlowList}

\title{Preprocess a list of flowSets to add empty flowFrames for patients 
       in which certain aliquots are missing }
\description{
The \code{preProcessFlowList} function takes in a list of flowSets as inputs
and adds empty flowFrames for patients in which data from  certain
aliquots are missing.  The function is to be used to modify datasets to
conform to a format that can be utilized by the multiple sample
\code{\link[flowQ:qaProcess-class]{qaProcess}} functions. 
}

\usage{
 preProcessFlowList(flowList)
}

\arguments{
  \item{flowList}{ A list of \code{flowSet} in which \code{flowFrames} for
  certain patients are missing}
}
 
\details{  
  This function can be used to modify a list of flowSets to fill in
  empty \code{flowFrames} for patients for which samples are missing.

}

\value{
  
  A \code{list} of objects of class \code{flowSet}.  
}

\author{ Nishant Gopalakrishnan }


\seealso{
  \code{\link[flowQ]{locateDuplicatedParameters}}
}




\keyword{ dynamic }



