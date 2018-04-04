gmmevents = function(time,identity,location,global_ids=NULL,verbose=TRUE, splitGroups=TRUE){
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
		
		DATA_LOC = DATA[location_indices,,drop=FALSE]
		
		igfdresult = infer_graph_from_datastream_mmVB(DATA_LOC,verbose)
		output=igfdresult[[1]]
		gmm=igfdresult[[2]]
		
		B_LOC = output$B_hard_incidence_matrix;
		gathering_events_LOC = dim(B_LOC)[1]

		
		B_aux = matrix(0,nrow=gathering_events_LOC,ncol=total_individuals_global)
		B_aux[,IM_current_DATA$global_indices[match(1:ncol(B_LOC),IM_current_DATA$current_indices)]] = B_LOC
		colnames(B_aux)=IM_current_DATA$original_indices

		
		if (splitGroups==TRUE) {
			T.t <- output$EVENT_TIMES
			T.t <- data.frame(Start=T.t[,1],End=T.t[,2],stringsAsFactors=FALSE)
			time.t <- DATA_LOC[,1]
			identity.t <- DATA_LOC[,2]
	
			while (TRUE) {
	
				T.t$gbi.id <- 1:nrow(T.t)
				T.t$Order <- rank(T.t$Start)
		
				T.t$Unique <- NA
				for (i in 1:(max(T.t$Order)-1)) {
					T.t$Unique[which(T.t$Order==i)] <- T.t$End[which(T.t$Order==i)] < min(T.t$Start[which(T.t$Order>i)])
				}
		
		
				event <- rep(NA,length(time.t))
				cur.event <- vector()
				for (i in 1:length(time.t)) {
					if (any(time.t[i] %in% T.t$Start)) {
						cur.event <- c(cur.event,T.t$Order[which(T.t$Start == time.t[i])])
					}
					event[i] <- cur.event[length(cur.event)]
					if (i < length(time.t)) {
						if (any(time.t[i] %in% T.t$End) & time.t[i] < time.t[i+1]) {
							cur.event <- cur.event[which(cur.event != T.t$Order[which(T.t$End == time.t[i])])]
						}
					} else {
						if (any(time.t[i] %in% T.t$End)) {
							cur.event <- cur.event[which(cur.event != T.t$Order[which(T.t$End == time.t[i])])]
						}
					}
				}

				if (any(T.t$Unique == FALSE, na.rm=TRUE)) {
					g <- min(T.t$Order[which(T.t$Unique==FALSE)])
					
					current.start <- T.t$Start[which(T.t$Order==g)]
					current.end <- 	T.t$End[which(T.t$Order==g)]
					next.start <- T.t$Start[which(T.t$Order == g+1)]
					new.end <- max(time.t[which(time.t >= current.start & time.t < next.start)])
					event[which(time.t > new.end & event == g)] <- NA
					new.start1 <- min(time.t[which(is.na(event))])
					new.end1 <- max(time.t[which(is.na(event))])
		
					T.t$Start[which(T.t$Order==g)] <- current.start
					T.t$End[which(T.t$Order==g)] <- new.end
					old.group <- T.t$gbi.id[which(T.t$Order==g)]
					T.t <- rbind(T.t, data.frame(Start=new.start1,End=new.end1,gbi.id=nrow(T.t)+1,Order=NA,Unique=NA))
					B_aux <- rbind(B_aux,rep(0,ncol(B_aux)))
					new.group <- nrow(B_aux)
				
					new.group.ids <- table(identity.t[which(is.na(event))])
					B_aux[new.group,] <- 0
					B_aux[new.group,IM_current_DATA$global_indices[which(IM_current_DATA$current_indices %in% as.numeric(names(new.group.ids)))]] <- new.group.ids
				
					old.group.ids <- table(identity.t[which(event==g)])
					B_aux[old.group,] <- 0
					B_aux[old.group,IM_current_DATA$global_indices[which(IM_current_DATA$current_indices %in% as.numeric(names(old.group.ids)))]] <- old.group.ids
				} else {
					break;
				}
			}
			output$EVENT_TIMES <- as.matrix(T.t[,c(1,2)])
		}
		
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

	
	
	
	
		