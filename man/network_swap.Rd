\name{network_swap}
\alias{network_swap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform one (or more) random swap
}
\description{
Performs one (or more) random swap on the data and re-calculates network, returning both the new network and the data stream
}
\usage{
network_swap(association_data, data_format = "GBI", swaps = 1, 
	association_index = "SRI", association_matrix = NULL, 
	identities = NULL, which_identities = NULL, times = NULL, locations = NULL, 
	which_locations = NULL, start_time = NULL, end_time = NULL, classes = NULL, 
	which_classes = NULL, days = NULL, within_day = FALSE, 
	within_location = FALSE, within_class = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{association_data}{
a \code{K x N} matrix of K groups (observations, gathering events, etc.) and N individuals (all individuals that are present in at least one group) OR a \code{K x N x N} array of sampling periods.
}
  \item{data_format}{
\code{"GBI"} expect a group by individual matrix, \code{"SP"} Expect a sampling periods array
}
  \item{swaps}{
number of swaps (default = 1000)
}
  \item{association_index}{
\code{"SRI"} Simple ratio index, \code{"HWI"} Half-weight index (more to come)
}
  \item{association_matrix}{
provide a starting association matrix (see details)
}
  \item{identities}{
N vector of identifiers for each individual (column) in the group by individual matrix
}
  \item{which_identities}{
vector of identities to include in the network (subset of identities)
}
  \item{times}{
K vector of times defining the middle of each group/event
}
  \item{locations}{
K vector of locations defining the location of each group/event
}
  \item{which_locations}{
vector of locations to include in the network (subset of locations)
}
  \item{start_time}{
element describing the starting time for inclusion in the network (useful for temporal analysis)
}
  \item{end_time}{
element describing the ending time for inclusion in the network (useful for temporal analysis)
}
  \item{classes}{
N vector of types or class of each individual (column) in the group by individual matrix (for subsetting)
}
  \item{which_classes}{
vector of class(es)/type(s) to include in the network (subset of classes)
}
  \item{days}{
K vector of day stamp for each event (can be integer or string representing any period of time)
}
  \item{within_day}{
if \code{TRUE} then permutations will be done within the time periods
}
  \item{within_location}{
if \code{TRUE} then permutations will be done within the given locations
}
  \item{within_class}{
if \code{TRUE} then permutations will be done within the given classes
}
}
\details{
Performs one or more permutation swaps on the group by individual matrix as given by Whitehead (2008).  In order to save on memory use, this function computers the number of swaps and returns the association matrix and the data stream resulting from these, thus not needing to create a large stack of networks to store each permutation. This can then be implemented in a loop as shown in the example below. Note that this method is quite a bit slower than the network_permutation function.

This implementation allows permutations (swaps) to be restricted to within any of three classes.  Though each class is labelled, the function is flexible.  Hence, days can represent any time period (months, hours, etc.). However, unlike the network_permutation, the subsetting of the data must be done outside of this function (for reasons that might be obvious) - see the example below.
}
\value{
Returns a list containing an \code{N x N} matrix with the dyadic association rates of each pair of individuals after performing the swaps, and the \code{N x N} data stream post-swap, as two list elements.
}
\references{
Whitehead (2008) \emph{Analyzing Animal Societies}
}
\author{
Damien R. Farine
}

\examples{

# load data
data("group_by_individual")
data("times")

# subset GBI (to reduce run time of the example)
gbi <- gbi[,1:80]

# calculate network for data based on morning associations
network <- get_network(gbi, association_index="SRI", 
	times=times, start_time=0, end_time=max(times)/2)

# perform 100 permutations and calculate the coefficient
# of variance after each permutation.
# note that the subsetting is done outside of the function
library(raster)
cvs <- rep(NA,100)
network_perm = list(network,gbi[which(times <= max(times)/2),])
hours <- floor(times/3600)[which(times <= max(times)/2)]
for (i in 1:100) {
	network_perm <- network_swap(network_perm[[2]], swaps=1, 
		association_matrix=network_perm[[1]], days=hours, 
		within_day=TRUE)
	cvs[i] <- cv(network_perm[[1]])
}

# plot the results with the original network as a red dot
plot(cvs,pch=20,cex=0.5)
points(0,cv(network),cex=1,pch=20,col="red")

}