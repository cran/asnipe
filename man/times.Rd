\name{times}
\alias{times}
\docType{data}
\title{
Observation Time for each Flock Contained in the Group by Individual data
}
\description{
The start time for each flock in the group by individual data since the first flock observed.
}
\usage{data("times")}
\format{
  An integer array with the time (in seconds since the first flock of the day).
}
\source{
Farine, D.R., Garroway, C.J., Sheldon, B.C. (2012) Social Network Analysis of mixed-species flocks: exploring the structure and evolution of interspecific social behaviour. Animal Behaviour 84: 1271-1277.
}
\examples{
data("times")
data("group_by_individual")
rownames(gbi) <- times
}
\keyword{datasets}
