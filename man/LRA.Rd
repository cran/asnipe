\name{LRA}
\alias{LRA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dyadic Lagged Association Rate
}
\description{
Calculate lagged association rate \code{g(tau)} from Whitehead (2008) for each dyad individually
}
\usage{
LRA(group_by_individual, times, timejump, output_style = 1, min_time = NULL, 
	max_time = NULL, identities = NULL, which_identities = NULL, locations = NULL, 
	which_locations = NULL, start_time = NULL, end_time = NULL, classes = NULL, 
	which_classes = NULL, association_rate = TRUE)
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
  \item{output_style}{
either 1 or 2, see details
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
  \item{association_rate}{
calculate lagged rate of association (see details)
}
}
\details{
Calculates the dyadic lagged association rate.  The lagged rate of association incorporates the number of observations of each individuals as a simple ratio index within each time period, leading to a better estimation of the assocation rate for data where many observations of individuals can be made within a single time period.
}
\value{
If \code{output_style == 1} then a stack of matrices is returned that is \code{N x N x tau}.
If \code{output_style == 2} then a dataframe is returned containing the focal ID, associate, \code{tau}, and lagged association rate.
}
\references{
Expanded from Whitehead (2008)
}
\author{
Damien R. Farine
}

\examples{

data("group_by_individual")
data("times")
data("individuals")

## calculate lagged association rate
lagged_rates <- LRA(gbi,times,3600, classes=inds$SPECIES, which_classes="GRETI", output_style=2)

## do something (run a model, plot a surface, etc..)
}
