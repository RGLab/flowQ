\name{qData}
\alias{qData}
\docType{data}
\title{Data for demonstrating QA processes that compares dye parameters
  from each patient across multiple aliquots}

\description{
  A \code{list} of 8 \code{flowSet}. Each flowSet has data from four
  patients and each 
element in the list corresponds to an aliquot.

}
\usage{data(qData)}
\format{
  The object contains a list of 8 objects of  class \code{flowSet}, each
  composed 
  of 4  \code{flowFrames}.
  Each \code{flowFrame} correponds to a sample from an aliquot.
}

\details{
  This qData dataset contains flouresence information for samples 
 from four patients. The sample from each patient was separated into eight 
 aliquots and each aliquot was stained using a different combination of dyes.
 For each patient, some dyes appear more than once in multiple alquots. These 
 along with the Forward and Side scatter information can be used for Quality 
 control.

}

\keyword{datasets}
