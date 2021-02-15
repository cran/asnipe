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
	identities = NULL, which_identities = NULL, times = NULL, occurrences = NULL,
	locations = NULL,  which_locations = NULL, start_time = NULL, 
	end_time = NULL, classes = NULL, which_classes = NULL, 
	days = NULL, within_day = FALSE, within_location = FALSE, within_class = FALSE,
	enter_time = NULL, exit_time = NULL, symmetric=TRUE, trialSwap=TRUE)
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
  \item{enter_time}{
N vector of times when each individual entered the population
}
  \item{exit_time}{
N vector of times when each individual departed the population 
}
  \item{symmetric}{
Boolean to ensure that permutations maintain symmetry within sampling periods if using \code{data_format="SP"}
}
  \item{trialSwap}{
Boolean to include trial swaps (if true, then every attempted permutation is returned)
}

}
\details{
Performs permutations on the group by individual matrix as given by Whitehead (2008).  In order to save computing, only the recently swapped individuals are recalculated, hence why the association matrix of the original data can be provided or is recalculated.

This implementation allows permutations (swaps) to be restricted to within any of three classes.  Though each class is labelled, the function is flexible.  Hence, days can represent any time period (months, hours, etc.).

Trial swaps are implemented following Miklos and Podani (2004). Every permutation, a candidate swap is selected. If the selected portions of the data satisfy the baseline rules (e.g. the checkerboard pattern), then either the selection is attempted again \code{trialSwap = FALSE} or not \code{trialSwap = TRUE}. This should be set to TRUE, but the option for FALSE is provided for legacy analyses.

See get_network function for additional details on each field.
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

### USING TIMES, ETC.

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



#### DOUBLE PERMUTATION EXAMPLE (see Farine & Carter 2021)

## Load data
data("group_by_individual")
data("times")

# subset GBI (to reduce run time of the example)
gbi <- gbi[,1:40]

# Specify metric
metric <- "DEGREE"

# calculate observed network
network <- get_network(gbi, data_format="GBI", 
	association_index="SRI", times=times)

# Calculate observed metric (degree)
degrees <- rowSums(network)

# Do randomisation (as above, permutations should be >=1000)
networks.perm <- network_permutation(gbi, data_format="GBI",
	association_matrix=network, times=times, permutations=10)

# Now calculate the same metric on all the random networks
degrees.rand <- apply(networks.perm,1,function(x) { rowSums(x)})

# Now substract each individual's median from the observed
degree.controlled <- degrees - apply(degrees.rand,1,median)

#### Now use degree.controlled for any later test. For example, to related against a trait:

# Make a trait
trait <- rnorm(length(degree.controlled))

# get the coefficient of this:
coef <- summary(lm(degree.controlled~trait))$coefficients[2,3]

# Compare this to a node permutation 
# (here just randomising the trait values)
# note this should be done >= 1000 times
n.node.perm <- 10
coefs.random <- rep(NA, n.node.perm)

for (i in 1:n.node.perm) {
	trait.random <- sample(trait)
	coefs.random[i] <- summary(lm(degree.controlled~trait.random))$coefficients[2,3]
}

# calculate P value (note this is only one sided)
P <- sum(coef <= coefs.random)/n.node.perm


}