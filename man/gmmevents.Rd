\name{gmmevents}
\alias{gmmevents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer gathering events
}
\description{
Infer gathering events (groups or flocks) from a temporal datastream of observations, such as PIT tag data.
}
\usage{
gmmevents(time, identity, location, global_ids=NULL, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
The timestamp for the observation. Must be a real number (i.e. not a date or time format). See details below.
}
  \item{identity}{
The identify of the individual in each observation (can be a number or a string, e.g. PIT tag code).
}
  \item{location}{
The location of the observation (can be a number or a string).
}
  \item{global_ids}{
A vector of all the IDs in the study, used if consistency needs to be maintained across datasets.
}
  \item{verbose}{
Whether to print out progress and information.
}
}
\details{
The gmmevents function has 3 primary inputs: time, identity, and location:

The time must be a number representing a real valued time stamp. This can be number of seconds since the start of the day, number of seconds since the start of the study, hour of the day, Julian date, etc. The time stamps should represent a meaningful scale given the group membership definition - for example if an edge is the propensity to observe an individual at the same location on the same day, then time stamps should be the day value. In the example below, the time stamps are in seconds, because flocks of birds visit feeders over a matter of minutes and the group definition is being in the same flock (and these occur over seconds to minutes). The input must be numeric whole numbers.

The identity is the unique identifier for each individual. This should be consistent across all of the data sets. In the example here, PIT tags are given, but in broader analyses, we would convert these to ring (band) numbers because individuals can have different PIT tag numbers in the course of the study but never change ring numbers. The function will accept any string or numeric input.

The location is where the observation took place. This should reflect meaningful observation locations for the study. The function will accept any string or numeric inputs.

If the analysis is being conducted as part of a broader analysis in the same populations, it can be useful to get the results in a consistent form each time. In that case, the global_ids variable can be used to maintain consistency each time an analysis is run, regardless of which individuals were identified in the current input data. That is, the group by individual (gbi) matrix will include a column for every individual provided in global_ids.

Further notes on usage:

The gmm_events functions requires a few careful considerations. First, the amount of memory used is the square of the amount of data - so having many observations in a given location can run out of memory. With 16gb of RAM, generally up to 10,000 observations per location (per day - see next point) seems to be a safe limit.

The input data provide for each location should take into account any artificial gaps in the observation stream. For example, if there are gaps in data collection at a given location, then the location information provided into the gmm_events function should be split into two 'locations' to represent each continuous set of observations. For example, in the PIT tag data set provided there are 8 days of sampling. Providing gmm_events with only the location data from the original data will cause the gap between days to override any gaps between groups (or flocks) within a given day. To overcome this, instead of providing the gmm_events function with just the location, it is important to provide a location by day variable. This variable is then returned in the metadata and the information extracted out again (using strsplit - see example below). 

}
\value{
Returns a list containing three items:

The first item is the group by individual matrix (gbi), which is a matrix where each row is a gathering event - or group - and each column is an individual. Cells in the matrix have a value of 1 if the individual was observed in that gathering event, and 0 if not. 

The second item is a matrix containing three columns: the start time, end time, and location of each detected event. The number of rows in this events matrix matches the number of rows in the gbi file, and the rows correspond to one another (thus, row 3 of the gbi has the start and end times, and location, of row 3 of the events matrix).

The third item has the same structure as the group by individual matrix, but instead of being binary (0 or 1), it contains the number of observation of individuals in each event.
}
\references{
Psorakis, I., Roberts, S. J., Rezek, I., & Sheldon, B. C. (2012). Inferring social network structure in ecological systems from spatio-temporal data streams. Journal of the Royal Society Interface, 9(76), 3055-3066. doi:10.1098/Rsif.2012.0223 \cr
Psorakis, I., Voelkl, B., Garroway, C. J., Radersma, R., Aplin, L. M., Crates, R. A., Culina, A., Farine, D. R., Firth, J.A., Hinde, C.A., Kidd, L.R., Milligan, N.D., Roberts, S.J., Verhelst, B., Sheldon, B. C. (2015). Inferring social structure from temporal data. Behavioral Ecology and Sociobiology, 69(5), 857-866. doi:10.1007/s00265-015-1906-0

}
\author{
Ioannis Psorakis (original code) \cr
Julian Evans (R implementation) \cr
Damien R. Farine (current implementation)
}
\examples{

\donttest{

	library(asnipe)
	data("identified_individuals")
	
	# Create unique locations in time
	identified_individuals$Loc_date <- 
		paste(identified_individuals$Location,
		identified_individuals$Date,sep="_")

	# Provide global identity list (including individuals 
	# not found in these data, but that need to be included).
	# Not including this will generate gbi with only the
	# individuals provided in the data set (in this case 151
	# individuals)
	global_ids <- levels(identified_individuals$ID)

	# Subset input data (Day 1)
	identified_individuals <- identified_individuals[which(identified_individuals$Date < 2),]
	
	# Generate GMM data
	gmm_data <- gmmevents(time=identified_individuals$Time,
		identity=identified_individuals$ID,
		location=identified_individuals$Loc_date,
		global_ids=global_ids)

	# Extract output
	gbi <- gmm_data$gbi
	events <- gmm_data$metadata

	# Can also subset gbi to only individuals observed
	# in the dataset to give same answer as if 
	# global_ids had not been provided
	gbi <- gbi[,which(colSums(gbi)>0)]

	# Split up location and date data
	tmp <- strsplit(events$Location,"_")
	tmp <- do.call("rbind",tmp)
	events$Location <- tmp[,1]
	events$Date <- tmp[,2]

}

}