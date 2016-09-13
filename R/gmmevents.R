gmmevents = function(time,identity,location,global_ids=NULL,verbose=TRUE){
	# Ioannis Psorakis, Stephen J. Roberts, Iead Rezek and Ben Sheldon
	# Inferring social network structure in ecological systems from spatio-temporal data streams
	# University of Oxford 2011
	# Translated to R by Julian Evans <jevansbio@gmail.com>
	#
	# Packaged to asnipe and maintained by Damien Farine <dfarine@orn.mpg.de>

	#require (MASS)
	#require (Matrix)

	oldw <- getOption("warn")
	options(warn = -1)
	
	if (length(time) != length(identity)) {
		stop("times and identities are not the same length")
	}

	if (length(location) != length(identity)) {
		stop("locations and identities are not the same length")
	}

	if (length(time) != length(location)) {
		stop("times and locations are not the same length")
	}

	if (is.null(global_ids)) {
		global_ids <- unique(identity)
		global_ids <- sort(global_ids)
		if (verbose) {
			print("Automatically generating and sorting IDs")
			print(sprintf("Detected %d individuals",length(global_ids)))
		}
	}

	DATA <- data.frame(time,identity,location,stringsAsFactors=FALSE)
	
	if(length(global_ids)>1){
		unique_indices = sort(unique(DATA[,2]))
		total_individuals_current_DATA = length(unique_indices)
		unique_indices_global = unique(global_ids)

		if (length(unique_indices_global) != length(global_ids)) {
			stop("Duplicate IDs provides in global_ids")
		}

		total_individuals_global = length(unique_indices_global)
		
		IM_current_DATA<-data.frame(unique_indices_global,c(1:total_individuals_global))
		names(IM_current_DATA)<-c("original_indices","global_indices")
		IM_current_DATA$current_indices=rep(NA,nrow(IM_current_DATA))
		IM_current_DATA$current_indices[IM_current_DATA$original_indices%in%unique_indices]<-c(1:total_individuals_current_DATA)
		
		locations_global <- sort(unique(location))
		LO_current_DATA<-data.frame(original.name=as.character(locations_global),location.index=c(1:length(locations_global)),stringsAsFactors=FALSE)
		
		DATA[,2] <- IM_current_DATA$current_indices[match(DATA[,2],IM_current_DATA[,1])]
		#DATA[,2]<-sapply(1:length(DATA[,2]),function (x) IM_current_DATA[IM_current_DATA[,1]==DATA[x,2],3])

		DATA[,3] <- LO_current_DATA$location.index[match(DATA[,3],LO_current_DATA[,1])]
		DATA<-as.matrix(DATA)

	} else {
		stop("Too few individuals provided")
	}
	
	unique.locations <- LO_current_DATA$location.index
	total_locations = length(unique.locations)
	
	B = vector()
	T = vector()
	L = vector()

	for (location_index in c(1: total_locations)){
		
		if (verbose) {
			print(sprintf("Starting analysis of location: %s", LO_current_DATA$original.name[location_index]))
		}

		location_indices = DATA[,3] == unique.locations[location_index]
		if (sum(location_indices) == 0){
			next()
		}
		
		if (length(unique(DATA[location_indices,2]))<2){
			next()
		}
		
		DATA_LOC = DATA[location_indices,]
		
		igfdresult = infer_graph_from_datastream_mmVB(DATA_LOC,verbose)
		output=igfdresult[[1]]
		gmm=igfdresult[[2]]
		
		B_LOC = output$B_hard_incidence_matrix;
		gathering_events_LOC = dim(B_LOC)[1]
		
		B_aux = matrix(0,nrow=gathering_events_LOC,ncol=total_individuals_global)
		B_aux[,IM_current_DATA$global_indices[match(1:ncol(B_LOC),IM_current_DATA$current_indices)]] = B_LOC
		colnames(B_aux)=IM_current_DATA$original_indices
		
		B <- rbind(B,B_aux)
		T <- rbind(T,output$EVENT_TIMES)
		L <- c(L,rep(LO_current_DATA$original.name[location_index],nrow(output$EVENT_TIMES)))

	}
	T <- data.frame(Start=T[,1],End=T[,2],Location=L,stringsAsFactors=FALSE)
	gbi <- B
	gbi[gbi>0] <- 1
	BT=list(gbi=gbi,metadata=T,B=B)
	options(warn = oldw)
	return(BT)
}

	
	
	
	
		