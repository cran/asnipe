network_permutation <- function(association_data, data_format = "GBI", permutations=1000, returns=1, association_index = "SRI", association_matrix = NULL, identities = NULL, which_identities = NULL, times = NULL, occurrences = NULL, locations = NULL, which_locations = NULL, start_time = NULL, end_time = NULL, classes = NULL, which_classes = NULL, days = NULL, within_day = FALSE, within_location = FALSE, within_class = FALSE) {
	
	#### CHECK INPUTS
	if (is.null(association_data)) { stop("No association_data data!") }
	if (length(dim(association_data)) != 2 & data_format=="GBI") { stop("Invalid dimensions for association_data") }
	if (length(dim(association_data)) != 3 & data_format=="SP") { stop("Invalid dimensions for association_data") }
	if ((length(identities) != ncol(association_data) & !is.null(identities)) == TRUE) { stop("Length of identities does not match number of individuals") }
	if ((length(times) != nrow(association_data) & !is.null(times)) == TRUE) { stop("Length of times does not match number of groups") }
	if ((length(occurrences[1,]) != nrow(association_data) & !is.null(occurrences)) == TRUE) { stop("Number of occurrence periods does not match number of sampling periods") }
	if ((length(occurrences[,1]) != ncol(association_data) & !is.null(occurrences)) == TRUE) { stop("Number of individuals in occurrences does not match number of individuals in sampling periods") }
	if ((length(locations) != nrow(association_data) & !is.null(locations)) == TRUE) { stop("Length of locations does not match number of groups") }
	if ((length(classes) != ncol(association_data) & !is.null(classes)) == TRUE) { stop("Length of classes does not match number of individuals") }
	if ((!is.null(which_identities) & is.null(identities)) == TRUE) { stop("Cannot apply which_identities without identities data") }
	if ((!is.null(which_locations) & is.null(locations)) == TRUE) { stop("Cannot apply which_locations without locations data") }
	if ((!is.null(start_time) & is.null(times)) == TRUE) { stop("Cannot apply start_time without times data") }
	if ((!is.null(end_time) & is.null(times)) == TRUE) { stop("Cannot apply end_time without times data") }
	if ((!is.null(which_classes) & is.null(classes)) == TRUE) { stop("Cannot apply which_class without classes data") }
	if (!any(association_index %in% c("SRI","HWI"))) { stop("Unknown association_index") }
	if (within_day & is.null(days)) { stop("Cannot constrict within days if days are not supplied") }
	if (within_location & is.null(locations)) { stop("Cannot constrict within location if locations are not supplied") }
	if (within_class & is.null(classes)) { stop("Cannot constrict within class if classes are not supplied") }
	if (!is.null(locations)) { locations <- as.matrix(locations) } # Fixes bug with data frames

	#### SUBSET THE DATA
	# By identity
	if (!is.null(which_identities)) {
		if (data_format=="GBI") association_data <- association_data[,which(identities %in% which_identities)]
		if (data_format=="SP") association_data <- association_data[,which(identities %in% which_identities),which(identities %in% which_identities)]
		identities <- identities[which(identities %in% which_identities)]
	}
	
	# By time
	if (!is.null(start_time) & is.null(end_time)) { end_time <- max(times) }
	if (!is.null(end_time) & is.null(start_time)) { start_time <- min(times) }
	if (!is.null(start_time) & !is.null(end_time)) {
		subs <- which(times >= start_time & times <= end_time)
		if (data_format=="GBI") association_data <- association_data[subs,]
		if (data_format=="SP") association_data <- association_data[subs,,]
		locations <- locations[subs]
		times <- times[subs]
	}
	
	# By location
	if (!is.null(which_locations)) {
		subs <- which(locations %in% which_locations)
		if (data_format=="GBI") association_data <- association_data[subs,]
		if (data_format=="SP") association_data <- association_data[subs,,]
		locations <- locations[subs]
		times <- times[subs]
	}
	
	# By class
	if (!is.null(which_classes)) {
		if (data_format=="GBI") association_data <-  association_data[,which(classes %in% which_classes)]
		if (data_format=="SP") association_data <- association_data[,which(classes %in% which_classes),which(classes %in% which_classes)]
		identities <- identities[which(classes %in% which_classes)]
	}
	
	if (!within_day) { days <- rep(1,nrow(association_data)) }
	if (!within_location) { locations <- rep(1,nrow(association_data)) }
	if (!within_class) { classes <- rep(1,ncol(association_data)) }
	
	
	#### GENERATE NETWORK IF REQUIRED
	### Calculate Network
	

		
	if (is.null(association_matrix)) {
	
		cat(paste("No association matrix provided, generating ", ncol(association_data), " x ", ncol(association_data), " matrix\n"))
		fradj_sorted <- get_network(association_data,data_format=data_format, association_index = "SRI", identities = identities, 
					which_identities = which_identities, times = times, locations = locations, 
					which_locations = which_locations, start_time = start_time, end_time = end_time, 
					classes = classes, which_classes = which_classes)

	} else {
		fradj_sorted <- association_matrix
	}

	if (!is.null(identities)) {
		colnames(fradj_sorted) <- identities
		rownames(fradj_sorted) <- identities
	}	


	#### DO PERMUTATIONS
	cat(paste("Starting permutations, generating ", ceiling(permutations/returns) , " x ", ncol(association_data), " x ", ncol(association_data) , " matrix"))
	association_data_perm <- association_data
	fradj_sorted2 <- fradj_sorted
	n_inds <- ncol(association_data_perm)
	fradj_sorted_perm <- array(0, c(ceiling(permutations/returns),n_inds,n_inds))
	
	# Calculate network
	do.SR_perm <- function(GroupBy,input){
			tmp <- input[ ,GroupBy] + input
			x <- colSums(tmp==2)
			yab <- colSums(tmp==1)
			if (association_index == "SRI") {
				out <-  x / (x + yab)
			} else if (association_index == "HWI") {
				out <- x / (x + 0.5*yab)
			}
		return(out)
	}

	do.SR_perm.times <- function(GroupBy,input,times) {
		tmp <- input[ ,GroupBy] + input
		x <- colSums(tmp==2)
		yab <- apply(tmp,2,function(x) { sum(table(times[x==1])==2) })
		y <- colSums(tmp==1)-(2*yab)
		if (association_index == "SRI") {
			out <- x / (x + y + yab)
		} else if (association_index == "HWI") {
			out <- x / (x + y + 0.5*yab)
		}
		return(out)
	}
	
	do.SR2_perm <- function(i, a, association_index) {
		# how many times 1 seen together with all others
		x <- apply(a[,i,],2,sum)

		# how many times 1 but not others in a sampling period and vice versa
		n <- apply(a,1,rowSums)
		n[n>0] <- 1
		seen <- t(apply(n,1,function(x) x-n[i,]))
		ya <- rowSums(seen<0)
		yb <- rowSums(seen>0)

		# how many times 1 and others seen but not together
		seen <- t(apply(n,1,function(x) x+n[i,]))
		yab <- rowSums(seen>1) - x
		
		if (association_index == "SRI") {
			out <- x / (x + ya + yb + yab)
		} else if (association_index == "HWI") {
			out <- x / (x + ya + yb + 0.5*yab)
		}
		return(out)
	}

	do.SR2_perm.occurrences <- function (i, a, association_index, occurrences) {
		# how many times 1 seen together with all others
		x <- apply(a[,i,],2,sum)

		# how many times 1 but not others in a sampling period and vice versa
		seen <- sweep(occurrences,2,occurrences[i,],"+")
		yab <- rowSums(seen==2)-x
		ya_b <- rowSums(seen==1)
				
		if (association_index == "SRI") {
			out <- x / (x + ya_b + yab)
		} else if (association_index == "HWI") {
			out <- x / (x + ya_b + 0.5*yab)
		}
		return(out)
	}

	count <- 1
	##  GET PERMUTATION MATRICES
	for (n in c(1:permutations)) {
		repeat {
			a <- which(association_data_perm>0,arr.ind=TRUE)
			s <- sample(1:nrow(a),1)
			first <- a[s,]
			if (data_format=="GBI") {
				a <- a[which(locations[as.numeric(a[,1])] == locations[as.numeric(a[s,1])] & days[as.numeric(a[,1])] == days[as.numeric(a[s,1])] & classes[as.numeric(a[,2])] == classes[as.numeric(a[s,2])]),,drop=FALSE]
				second <- a[sample(1:nrow(a),1),]
				if (first[1]!=second[1]&first[2]!=second[2]&(sum(association_data_perm[first[1],first[2]]) > 0 & sum(association_data_perm[second[1],second[2]]) > 0) & (sum(association_data_perm[second[1],first[2]]) == 0 & sum(association_data_perm[first[1],second[2]]) == 0)) { break; }
			}
			if (data_format=="SP") {
				a <- a[which(a[,1] == a[s,1] & classes[as.numeric(a[,2])] == classes[as.numeric(a[s,2])] & classes[as.numeric(a[,3])] == classes[as.numeric(a[s,3])]),,drop=FALSE]
				second <- a[sample(1:nrow(a),1),]
				if (first[2]!=second[2]&first[3]!=second[3]&(sum(association_data_perm[first[1],first[2],first[3]]) > 0 & sum(association_data_perm[second[1],second[2],second[3]]) > 0) & (sum(association_data_perm[second[1],second[2],first[3]]) == 0 & sum(association_data_perm[first[1],first[2],second[3]]) == 0)) { break; }
			}
				
		}
		
		if (data_format=="GBI") {
			association_data_perm[second[1],first[2]] <- association_data_perm[first[1],first[2]]
			association_data_perm[first[1],second[2]] <- association_data_perm[second[1],second[2]]
			association_data_perm[first[1],first[2]] <- 0
			association_data_perm[second[1],second[2]] <- 0
		
			if (data_format=="GBI" & is.null(times)) {
				tmp1 <- do.SR_perm(first[2],association_data_perm)
				tmp2 <- do.SR_perm(second[2],association_data_perm)
			} else {
				tmp1 <- do.SR_perm.times(first[2],association_data_perm,times)
				tmp2 <- do.SR_perm.times(second[2],association_data_perm,times)
			}
		}
		if (data_format=="SP") {
			association_data_perm[first[1],second[2],first[3]] <- association_data_perm[first[1],first[2],first[3]]
			association_data_perm[first[1],first[2],second[3]] <- association_data_perm[first[1],second[2],second[3]]
			association_data_perm[first[1],first[2],first[3]] <- 0
			association_data_perm[first[1],second[2],second[3]] <- 0
		
			if (is.null(occurrences)) {
				tmp1 <- do.SR2_perm(first[2],association_data_perm, association_index)
				tmp2 <- do.SR2_perm(second[2],association_data_perm, association_index)
			} else {
				tmp1 <- do.SR2_perm.occurrences(first[2],association_data_perm, association_index, occurrences)
				tmp2 <- do.SR2_perm.occurrences(second[2],association_data_perm, association_index, occurrences)
			}
		}

		fradj_sorted2[,first[2]] <- tmp1
		fradj_sorted2[first[2],] <- tmp1
		fradj_sorted2[,second[2]] <- tmp2
		fradj_sorted2[second[2],] <- tmp2
		
		diag(fradj_sorted2) <- 0
		
		if ((n %% returns) == 0) {
			fradj_sorted_perm[count,,] <- fradj_sorted2
			count <- count + 1
		}
	}
	
	if ((n %% returns) != 0) {
		fradj_sorted_perm[count,,] <- fradj_sorted2
		count <- count + 1
	}
	return(fradj_sorted_perm)
}
