\name{get_sampling_periods}
\alias{get_sampling_periods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert group or individual data into sampling periods
}
\description{
Converts several different types of data storage into sampling periods for calculating or permuting networks
}
\usage{
get_sampling_periods(association_data, association_times, sampling_period, 
identities = NULL, location = NULL, within_locations = FALSE, 
data_format = c("gbi", "groups", "individuals"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{association_data}{
Can be either a group by individual matrix, a list containing group members in each element, or a two-column data frame with individual ID in the first column and group ID in the second column
}
  \item{association_times}{
Because sampling periods are inferred over time, each group must contain some time data (can be in any format, such as seconds, days, etc.). One time must be provided for each row of the association data.
}
  \item{sampling_period}{
The number of time periods over which data are combined (for example 10 days, 3600 seconds)
}
  \item{identities}{
Optional identities for each individual in the dataset
}
  \item{location}{
If spatial disaggregation need to be maintained, samping periods can be calculated per time per location
}
  \item{within_locations}{
Flag whether to include location information
}
  \item{data_format}{
Format of the input data
}
}
\details{
This function will calculate an association matrix for each sampling period. If locations are included, these will be treated independently.
}
\value{
Returns a \code{K x N x N} stack of matrices, where each \code{N x N} slice is an association matrix. Row names and Column names of these slices are given the identity where available. The K slice names are given either the time or time_location for each sampling period.
}
\author{
Damien R. Farine
}
\examples{
## define group memberships (these would be read from a file)
individuals <- data.frame(ID=c("C695905","H300253","H300253",
	"H300283","H839876","F464557","H300296","H300253",
	"F464557","H300296","C695905","H300283","H839876"),
	GROUP=c(1,1,2,2,2,3,3,4,5,5,6,6,6))

## create a time column
individuals <- cbind(individuals,
	DAY=c(1,1,1,1,1,2,2,2,3,3,3,3,3))

SPs <- get_sampling_periods(individuals[,c(1,2)],
	individuals[,3],1,data_format="individuals")
	
## define group memberships (these would be read from a file)
groups <- list(G1=c("C695905","H300253"),
	G2=c("H300253","H300283","H839876"),
	G3=c("F464557","H300296"),
	G4=c("H300253"),
	G5=c("F464557","H300296"),
	G6=c("C695905","H300283","H839876"))

## create a time variable
days <- c(1,1,2,2,3,3)

SPs <- get_sampling_periods(groups,
	days,1,data_format="groups")

}
