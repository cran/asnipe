\name{network_permutation}
\alias{network_permutation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Permutation
}
\description{
Performs permutations on the data and calculates network for each step
}
\usage{
network_permutation(association_data, data_format = "GBI", permutations = 1000, 
	returns=1, association_index = "SRI", association_matrix = NULL, 
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
  \item{permutations}{
number of permutations (default = 1000)
}
  \item{returns}{
number of swaps to perform between each association matrix that is returned (default = 1)
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
Performs permutations on the group by individual matrix as given by Whitehead (2008).  In order to save computing, only the recently swapped individuals are recalculated, hence why the association matrix of the original data can be provided or is recalculated.

This implementation allows permutations (swaps) to be restricted to within any of three classes.  Though each class is labelled, the function is flexible.  Hence, days can represent any time period (months, hours, etc.).
}
\value{
Returns a \code{p x N x N} stack of matrices with the dyadic association rates of each pair of individuals after each swap or after a number of swaps, where p = ceiling(permutations/returns)
}
\references{
Whitehead (2008) \emph{Analyzing Animal Societies}
}
\author{
Damien R. Farine
}

\examples{

data("group_by_individual")
data("times")

# subset GBI (to reduce run time of the example)
gbi <- gbi[,1:80]

## define to 2 x N x N network to hold two association matrices
networks <- array(0, c(2, ncol(gbi), ncol(gbi)))

## calculate network for first half of the time
networks[1,,] <- get_network(gbi, data_format="GBI", 
	association_index="SRI", times=times, start_time=0, 
	end_time=max(times)/2)
networks[2,,] <- get_network(gbi, data_format="GBI",
	association_index="SRI", times=times, 
	start_time=max(times)/2, end_time=max(times))

## calculate the weighted degree
library(sna)
deg_weighted <- degree(networks,gmode="graph", g=c(1,2), ignore.eval=FALSE)

## perform the permutations constricting within hour of observation
## note permutations are limited to 10 to reduce runtime
network1_perm <- network_permutation(gbi, data_format="GBI",
	association_matrix=networks[1,,], times=times, start_time=0, 
	end_time=max(times)/2, days=floor(times/3600), within_day=TRUE,
	permutations=10)
network2_perm <- network_permutation(gbi, data_format="GBI",
	association_matrix=networks[2,,], times=times, 
	start_time=max(times)/2, end_time=max(times), days=floor(times/3600), within_day=TRUE,
	permutations=10)

## calculate the weighted degree for each permutation
deg_weighted_perm1 <- degree(network1_perm,gmode="graph", g=c(1:10), ignore.eval=FALSE)
deg_weighted_perm2 <- degree(network2_perm,gmode="graph", g=c(1:10), ignore.eval=FALSE)
detach(package:sna)

## plot the distribution of permutations with the original data overlaid
par(mfrow=c(1,2))
hist(colMeans(deg_weighted_perm1),breaks=100, 
	main=paste("P = ", 
	sum(mean(deg_weighted[,1]) < colMeans(deg_weighted_perm1))/ncol(deg_weighted_perm1)), 
	xlab="Weighted degree", ylab="Probability")
abline(v=mean(deg_weighted[,1]), col='red')
hist(colMeans(deg_weighted_perm2),breaks=100, 
	main=paste("P = ",
	sum(mean(deg_weighted[,2]) < colMeans(deg_weighted_perm2))/ncol(deg_weighted_perm2)), 
	xlab="Weighted degree", ylab="Probability")
abline(v=mean(deg_weighted[,2]), col='red')

}