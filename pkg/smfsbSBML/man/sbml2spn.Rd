\name{sbml2spn}
\alias{sbml2spn}
\title{Read an SBML file and parse into a smfsb SPN object}
\description{
This function reads and parses an SBML file into an SPN object suitable for simulation using the "smfsb" R package.
}
\usage{
sbml2spn(filename, verb=FALSE)
}
\arguments{
  \item{filename}{The path to the SBML file.}
  \item{verb}{If TRUE, print some basic parsing info to the console.}
}
\value{
An SPN object representing the a reaction network model.}
\seealso{\code{\link[smfsb]{StepGillespie}}, \code{\link[smfsb]{simTs}}}

\keyword{smfsb}
\keyword{libSBML}
\keyword{smfsbSBML}

