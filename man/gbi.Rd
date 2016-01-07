\name{gbi}
\alias{gbi}
\docType{data}
\title{
Detections of Individuals Forming Flocks at Bird Feeders
}
\description{
Dataset consisting of 151 individuals of 5 passerine species in Wytham Woods, UK: 78 blue tits (Cyanistes caeruleus), 7 coal tits (Periparus ater), 51 great tits (Parus major), 11 marsh tits (Poecile palustris), 3 nuthatches (Sitta europaea) and 1 individual of unknown species. Individuals were all fitted with individually-encoded passive integrated transponder (PIT) tags  that were logged by radio frequency identification (RFID) antennae fitted to each hole on regular sunflower feeders (we used unhusked sunflower seed). Data were collected from 4 feeders spaced approximately 300m over the course of one day. Feeders logged the presence of individuals at a sub-second resolution, and detections were assigned to flocks using a machine learning algorithm (a Gaussian Mixture Model).
}
\usage{data("group_by_individual")}
\format{
  Data are formatted in a group by individual matrix. Each row represents one flock, each column represents one individual.
}
\source{
Farine, D.R., Garroway, C.J., Sheldon, B.C. (2012) Social Network Analysis of mixed-species flocks: exploring the structure and evolution of interspecific social behaviour. Animal Behaviour 84: 1271-1277.
}
\examples{
data("group_by_individual")
str(gbi) # see the structure of the data
}
\keyword{datasets}
