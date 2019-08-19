get_sampling_periods <- function(association_data, association_times, sampling_period, identities = NULL, location = NULL, within_locations = FALSE, data_format = c("gbi","groups", "individuals"), return="SP") {


	if (is.numeric(association_times)) {
		period <- ceiling(association_times/sampling_period)
	} else {
		if (sampling_period == 1) {
			period <- as.character(association_times)
		} else {
			stop("Cannot calculate periods with sampling_period != 1 and non-numeric association_times")
		}
	}
	periods <- unique(period)
	
	if (within_locations == TRUE) {
		if (!is.null(location)) {
			locations <- unique(location)
			period <- paste(period,location,sep="_")
			periods <- unique(period)
		} else {
			stop("No location information provided")
		}
	}
	
	n_inds <- length(identities)
	if (is.null(identities)) {
		if (data_format == "groups") {
			n_inds <- length(unique(unlist(association_data, use.names = FALSE)))
			ids <- unique(unlist(association_data, use.names = FALSE))
		}
		if (data_format == "gbi") {
			n_inds <- ncol(association_data)
			if (!is.null(colnames(association_data))) { ids <- colnames(association_data);
			} else { ids <- c(1:ncol(association_data)); }
			association_data <- apply(association_data,1,function(x) ids[x>0])
		}
		if (data_format == "individuals") {
			n_inds <- length(unique(association_data[,1]))
			ids <- unique(association_data[,1])
		}
	}
	if (n_inds==0) stop("Error calculating number of individuals")
	
	if (return=="SP") {
	sampling_periods <- array(0,c(length(periods),n_inds,n_inds), dimnames = list(periods, ids, ids))
		
	for (i in periods) {
		
		# GBI format
		#if (data_format == "gbi") {
		#	if (sum(period==i)>1) {
		#		inds <- apply(association_data[period==i,],1,function(x) ids[which(x > 0)], as.list)
		#	} else {
		#		inds <- list(ids[which(association_data[period==i,] > 0)])
		#	}
		#}
		
		# individuals format
		if (data_format == "individuals") {
			inds_tmp <- association_data[period==i,]
			grps <- unique(inds_tmp[,2])
			inds <- list()
			for (j in 1:length(grps)) {
				inds[[j]] <- inds_tmp[which(inds_tmp[,2] %in% grps[j]),1]
			}
		}
		
		if (data_format == "groups" | data_format == "gbi") inds <- association_data[period==i]
		
		# inds is a list of associations
		for (j in 1:length(inds)) {
			sampling_periods[which(periods==i),which(ids %in% inds[[j]]),which(ids %in% inds[[j]])] <- 1
		}
		diag(sampling_periods[which(periods==i),,]) <- 0
	
	}
	
	return(sampling_periods)
	} else {
		if (sampling_period == 1 & within_locations == FALSE) {
		if (data_format == "individuals") {
			tab=table(association_data[,1],sampling_period=association_times)
			return(tab[,periods])
		}
		if (data_format == "groups") {
			tab <- matrix(0, n_inds, length(periods), dimnames=list(ids, periods))
			for (i in 1:length(periods)) {
				tab[which(ids %in% unlist(association_data[period==periods[i]])),i] <- 1
			}
			return(tab)
		}
		if (data_format == "gbi") {
			tab <- matrix(0, n_inds, length(periods), dimnames=list(ids, periods))
			for (i in 1:length(periods)) {
				tab[,i] <- tab[,i] + (colSums(association_data[period==periods[i],,drop=FALSE]) > 0)
			}
			return(tab)
		}
		} else {
			stop("Sorry not implemented for multiple locations or summarise sampling periods")
		}
	}
	
}
