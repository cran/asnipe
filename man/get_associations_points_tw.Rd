\name{get_associations_points_tw}
\alias{get_associations_points_tw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Group Membership using Time Window (please read warnings before using this method)
}
\description{
A time window approach to calculate group co-memberships.
}
\usage{
get_associations_points_tw(point_data, time_window = 180, which_days = NULL, 
	which_locations = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{point_data}{
dataframe of four columns: \code{Date  Time  ID  Location}.  This requirement is strict (see details).
}
  \item{time_window}{
window around each individual for calculating associations
}
  \item{which_days}{
subset of \code{Date} to include
}
  \item{which_locations}{
subset of \code{Locations} to include
}
}
\details{
Calculates an ego-centric approach to group membership (see warning).  For each detection, a group is created with and all associates within the time window at the same location are included.  

Input data must be of the following form:  Date is an integer for day (usually starting at 1 on the first day).  Time are the number of seconds elapsed from the start (continuous across all dates).  ID is a unique character string for each individual.  Location is a unique character string for each location.
}
\value{
Returns a list with three objects:
1. group by individual matrix (K rows by N columns)
2. an vector of times for each group
3. a vector of dates for each group
4. a vector of locations for each group
}
\author{
Damien R. Farine
}
\section{Warning }{
This method biases associations of dyads occuring in large groups because it creates one row in the association matrix for each detection of an individual. For this reason, this function should not be used (see also Psorakis et al. 2015 Behavioural Ecology & Sociobiology). One way to circumvent this is by including only groups centered around the focal individual when calculating associations. However, none of the functions in this package are implement this way.
}



\examples{

data("identified_individuals")

## calculate group_by_individual for first day at one location
group_by_individual <- get_associations_points_tw(identified_individuals, time_window=180, 
	which_days=1,which_locations="1B")

## split the resulting list
times <- group_by_individual[[2]]
dates <- group_by_individual[[3]]
locations <- group_by_individual[[4]]
group_by_individual <- group_by_individual[[1]]

}