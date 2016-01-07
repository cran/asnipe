\name{get_group_by_individual}
\alias{get_group_by_individual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert group or individual data into a group by individual matrix
}
\description{
Converts several different types of data storage into a group by individual matrix for calculating or permuting networks
}
\usage{
get_group_by_individual(association_data, identities = NULL, 
location = NULL, data_format = c("groups", "individuals"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{association_data}{
Can be either a group by individual matrix or a list containing group members in each element
}
  \item{identities}{
Optional identities for each individual in the dataset
}
  \item{location}{
Returns these spatial locations for each group
}
  \item{data_format}{
Format of the input data
}
}
\details{
This function will calculate an \code{K x N} matrix representing \code{K} groups and \code{N} individuals. If locations are included, these will be returned in the row names.
}
\value{
Returns a \code{K x N} matrix, where each \code{K} row is an group defined from the input data. Column names of the matrix are given the identity where available. The K row names are given either the time or time_location for each group.
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

gbi <- get_group_by_individual(individuals, 
	data_format="individuals")
	
## define group memberships (these would be read from a file)
groups <- list(G1=c("C695905","H300253"),
	G2=c("H300253","H300283","H839876"),
	G3=c("F464557","H300296"),
	G4=c("H300253"),
	G5=c("F464557","H300296"),
	G6=c("C695905","H300283","H839876"))

## create a time variable
days <- c(1,1,2,2,3,3)

gbi <- get_group_by_individual(groups, 
	data_format="groups")

}
