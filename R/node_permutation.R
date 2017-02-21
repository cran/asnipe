node_permutation <- function(network, permutations=1000, identities = NULL, which_identities = NULL, locations = NULL, which_locations = NULL, classes = NULL, which_classes = NULL, within_location = FALSE, within_class = FALSE) {
	
	#### CHECK INPUTS
	if (is.null(network)) { stop("No network data!") }
	if (length(dim(network)) != 2) { stop("Invalid dimensions for network data") }
	if (dim(network)[1] != dim(network)[2]) { stop("Invalid dimensions for network data") }
	if ((length(identities) != ncol(network) & !is.null(identities)) == TRUE) { stop("Length of identities does not match number of individuals") }
	if ((length(locations) != nrow(network) & !is.null(locations)) == TRUE) { stop("Length of locations does not match number of individuals") }
	if ((length(classes) != ncol(network) & !is.null(classes)) == TRUE) { stop("Length of classes does not match number of individuals") }
	if ((!is.null(which_identities) & is.null(identities)) == TRUE) { stop("Cannot apply which_identities without identities data") }
	if ((!is.null(which_locations) & is.null(locations)) == TRUE) { stop("Cannot apply which_locations without locations data") }
	if ((!is.null(which_classes) & is.null(classes)) == TRUE) { stop("Cannot apply which_class without classes data") }
	if (within_location & is.null(locations)) { stop("Cannot constrict within location if locations are not supplied") }
	if (within_class & is.null(classes)) { stop("Cannot constrict within class if classes are not supplied") }
	if (!is.null(locations)) { locations <- as.matrix(locations) } # Fixes bug with data frames


	# ADD DATA
	if (is.null(locations)) { locations <- rep(1,dim(network)[1]) }
	if (is.null(classes)) { classes <- rep(1,dim(network)[1]) }

	#### SUBSET THE DATA
	# By identity
	if (!is.null(which_identities)) {
		network <- network[which(identities %in% which_identities),which(identities %in% which_identities)]
		identities <- identities[which(identities %in% which_identities)]
	}
	
	# By location
	if (!is.null(which_locations)) {
		network <- network[which(locations %in% which_locations),which(locations %in% which_locations)]
		locations <- locations[which(locations %in% which_locations)]
	}
	
	# By class
	if (!is.null(which_classes)) {
		network <- network[which(classes %in% which_classes),which(classes %in% which_classes)]
		identities <- identities[which(classes %in% which_classes)]
	}
	
	if (!within_location) { locations <- rep(1,nrow(network)) }
	if (!within_class) { classes <- rep(1,ncol(network)) }
	
	
	if (!is.null(identities)) {
		colnames(fradj_sorted) <- identities
		rownames(fradj_sorted) <- identities
	}	


	#### DO PERMUTATIONS
	cat(paste("Starting permutations, generating ", permutations , " x ", ncol(network), " x ", ncol(network) , " matrix"))

	network.random <- array(NA, c(permutations,nrow(network),ncol(network)))

	classes.locations <- paste(classes,locations,sep="_")
	unique.classes.locations <- unique(classes.locations)

	for (i in 1:permutations) {
		
		for (j in 1:length(unique.classes.locations)) {
			
			sub <- which(classes.locations == unique.classes.locations[j])
			sub.net <- network[sub,sub,drop=FALSE]
			new.order <- sample(1:nrow(sub.net))
			sub.net.rand <- sub.net[new.order,new.order]
			network.random[i,sub,sub] <- sub.net.rand
		
		}

	}

	return(network.random)
}
