\name{txtFormatQAObject}
\alias{txtFormatQAObject}

\title{Formats results contained in QA process objects into a data frame with 
		appropriate row and column names}
\description{
Formats information contained in QA process objects into a data frame
which can then be easly written to a tab delimited file. Information
contained includes the values used in the outlier detection, information
whether the particular parameter passed the outlier  test and the
summary aggregator information. 
}
\usage{
  txtFormatQAObject(qp)
}

\arguments{
  \item{qp}{ A \code{qaProcess} object }
  
}
 


\value{
  
  A \code{data.frame} containing the values used in outlier detection,
  thresholds used, information regarding whether the parameter
  passed/failed the test as  well as the summary aggregator information
}

\author{ Nishant Gopalakrishnan }


\seealso{
  \code{\link[flowQ]{writeQATextReport}}
}


\examples{
  data(qData)
  p1 <- qaProcess.marginevents(qData[[1]], channels=c("FSC-A", "SSC-A"),
           outdir=tempdir(), name=paste("Panel", 1, sep="_"))
  res <- txtFormatQAObject(p1)
  res

}

\keyword{ dynamic }


 
