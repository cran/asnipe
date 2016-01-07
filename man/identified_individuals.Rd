\name{identified_individuals}
\alias{identified_individuals}
\docType{data}
\title{
Raw Observation Data of Individual Birds Feeding at Flocks
}
\description{
Contains the raw observation data, of which the first day was used to form the group by individual file. IDs correspond to TAG in the "individuals" data (note that some tags are error codes, which have not been removed, and thus do not occur in the individuals data).
}
\usage{data("identified_individuals")}
\format{
  Data frame containing 4 columns:
  	Date - The observation day (1 to 8, where days 1-2 are the first weekend, 3-4 the second weekend, etc..)
  	Time - The time in seconds since the very first observation
  	ID - The PIT tag code of the individual
  	Location - The location where the detection was made (1B, 1C, 1D, 1E)
}
\source{
Farine, D.R., Garroway, C.J., Sheldon, B.C. (2012) Social Network Analysis of mixed-species flocks: exploring the structure and evolution of interspecific social behaviour. Animal Behaviour 84: 1271-1277.
}
\examples{
data("identified_individuals")
head(identified_individuals)
table(identified_individuals$Location)
}
\keyword{datasets}
