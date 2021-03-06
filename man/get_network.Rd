\name{get_network}
\alias{get_network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating Weighted Network
}
\description{
Calculate a network from a group by individual matrix.  This function allows various levels of subsetting.
}
\usage{
get_network(association_data, data_format = "GBI", 
	association_index = "SRI", identities = NULL, 
	which_identities = NULL, times = NULL, occurrences = NULL,
	locations = NULL, which_locations = NULL, start_time = NULL, 
	end_time = NULL, classes = NULL, which_classes = NULL,
	enter_time = NULL, exit_time = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{association_data}{
a \code{K x N} matrix of K groups (observations, gathering events, etc.) and N individuals (all individuals that are present in at least one group) OR a \code{K x N x N} array of sampling periods.
}
  \item{data_format}{
\code{"GBI"} expect a group by individual matrix, \code{"SP"} Expect a sampling periods array
}
  \item{association_index}{
\code{"SRI"} Simple ratio index, \code{"HWI"} Half-weight index (more to come)
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
  \item{occurrences}{
N x S matrix with the occurrence of each individual in each sampling period (see details) containing only 0s and 1s
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
  \item{enter_time}{
N vector of times when each individual entered the population
}
  \item{exit_time}{
N vector of times when each individual departed the population 
}
}
\details{
Provides the ability to generate networks from one group by individual matrix and subsetting within the function.  This is particularly useful for generating several networks with different characteristics from the same group by individual matrix (for example networks from a given location or set of locations, or of a particular sex).

Including occurrence data is recommended when using sampling periods (not required for GBI data). If an individual is only observed alone in a sampling period, then it will not be included in the sampling period matrices (as these record only associations or interactions, not presence). Thus, a matrix containing N (for number of individuals) rows and S (for number of sampling periods) is required. See the get_sampling_periods function for help generating this matrix.

In some situations it is useful to calculate the network based only on the period in which each dyad overlapped within the population. In such cases, the entry_time and/or the exit_time variables can be given. These must given in the same format as the times variable, and all need to be in a format capable of doing time or date comparisons using > and < operators. The easiest is YYYYMMDD, whereas MMDDYYYY or DDMMYYYY will not work properly.

}
\value{
\code{N x N} matrix of association weights for each dyad.
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

## test if one predicts the other via a mantel test (must be loaded externally)
library(ape)
mantel.test(networks[1,,],networks[2,,])

## convert to igraph network and calculate degree of the first network
\dontrun{
library(igraph)
net <- graph.adjacency(networks[1,,], mode="undirected", diag=FALSE, weighted=TRUE)
deg_weighted <- graph.strength(net)
detach(package:igraph)


## alternatively package SNA can use matrix stacks directly
library(sna)
deg_weighted <- degree(networks,gmode="graph", g=c(1,2), ignore.eval=FALSE)
detach(package:sna)
}
}