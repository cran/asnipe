test_network_hypothesis <- function (association_data, type="GBI", randomisation="D-S", test=mean, random.networks=NULL, ..., plot=FALSE) {

	arguments <- list(...)

	if (type %in% c("GBI","SP")) {
		arguments.to.pass <- arguments[which(names(arguments) %in% names(formals(get_network)))]
		network <- do.call(get_network, c(list(association_data=association_data, data_format=type), arguments.to.pass))
	} else if (type == "Network") {
		network <- association_data
	} else {
		stop ("Unknown association_data type")
	}

	if (is.null(random.networks)) {
		
		if (randomisation=="D-S" & type %in% c("GBI","SP")) {

			arguments.to.pass <- arguments[which(names(arguments) %in% names(formals(network_permutation)))]
			random.networks <- do.call(network_permutation, c(list(association_data=association_data, association_matrix=network), arguments.to.pass))

		} else if (randomisation=="Node" & type == "Network") {

			arguments.to.pass <- arguments[which(names(arguments) %in% names(formals(node_permutation)))]
			random.networks <- do.call(node_permutation,c(list(network=association_data), arguments.to.pass))

		} else {
		
			stop("Unknown randomisation specified, or mismatch between association data type and randomisation type")

		}

	}

	
	observed <- test(network)
	random <- rep(NA, dim(random.networks)[1])

	for (i in 1:dim(random.networks)[1]) {

		random[i] <- test(random.networks[i,,])

	}

	if (plot==TRUE) {
	
		hist(random,breaks=50,col="lightgrey", main=paste("P =",round(sum(observed>=random)/length(random),3)))
		abline(v=observed, col="red")

	}

	return(list(observed=observed, random=random, P=sum(observed <= random)/length(random)))

}

