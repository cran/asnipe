\name{inds}
\alias{inds}
\docType{data}
\title{
Data on the Individual Birds Contained in the Group by Individual data
}
\description{
Information about the PIT tag number, ring number, species, and sex (where available) for each individual in the group by individual data. Each row represents one column in the group by individual file, and the order is maintained.
}
\usage{data("individuals")}
\format{
  Data frame containing:
  	TAG - A 10 character hexadecimal code unique to each individual
  	RING.NUMBER - A 7 character unique ring (or band) number for each individual
  	SPECIES - Each species, where BLUTI=blue tit, COATI=coal tit, GRETI=great tit, MARTI=marsh tit, and NUTHA=nuthatch
}
\source{
Farine, D.R., Garroway, C.J., Sheldon, B.C. (2012) Social Network Analysis of mixed-species flocks: exploring the structure and evolution of interspecific social behaviour. Animal Behaviour 84: 1271-1277.
}
\examples{
data("individuals")
data("group_by_individual")
colnames(gbi) <- inds$RING.NUMBER
}
\keyword{datasets}
