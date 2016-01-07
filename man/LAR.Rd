\name{LAR}
\alias{LAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Lagged Association Rate
}
\description{
Calculate lagged association rate \code{g(tau)} from Whitehead (2008)
}
\usage{
LAR(group_by_individual, times, timejump, min_time = NULL, max_time = NULL, 
	identities = NULL, which_identities = NULL, locations = NULL, 
	which_locations = NULL, start_time = NULL, end_time = NULL, classes = NULL, 
	which_classes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group_by_individual}{
a \code{K x N} matrix of K groups (observations, gathering events, etc.) and N individuals (all individuals that are present in at least one group)
}
  \item{times}{
K vector of times defining the middle of each group/event
}
  \item{timejump}{
step length for \code{tau}
}
  \item{min_time}{
minimum/starting value of \code{tau}
}
  \item{max_time}{
maximum/ending value of \code{tau}
}
  \item{identities}{
N vector of identifiers for each individual (column) in the group by individual matrix
}
  \item{which_identities}{
vector of identities to include in the network (subset of identities)
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
}
\details{
Calculate the lagged association rate for given timesteps.
}
\value{
Returns a matrix with Log(time) in the first column and the lagged association rate in the second
}
\references{
Whitehead (2008) \emph{Analyzing Animal Societies} section 5.5.1
}
\author{
Damien R. Farine
}

\examples{

data("group_by_individual")
data("times")
data("individuals")

## calculate lagged association rate for great tits
lagged_rates <- LAR(gbi,times,3600, classes=inds$SPECIES, which_classes="GRETI")

## plot the results
plot(lagged_rates, type='l', axes=FALSE, xlab="Time (hours)", ylab="LAR", ylim=c(0,1))
axis(2)
axis(1, at=lagged_rates[,1], labels=c(1:nrow(lagged_rates)))

}
