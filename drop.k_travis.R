drop.k <- function(runjags.object, dropvars, k=1, simulations=NA, ...){
	
	runjags.object <- runjags:::checkvalidrunjagsobject(runjags.object)
	if(runjags.object$sample==0)
		stop('Cannot perform a drop-k study on a simulation with no samples')
	
	if(missing(dropvars))
		stop('The name of the variable to be dropped must be specified as dropvars')
	
	# Some options that aren't allowed:
	if(any(c('model','datafunction','inits', 'monitor', 'data', 'n.chains', 'combine', 'datalist', 'initlist', 'export.cluster', 'test') %in% names(list(...)))){
		stop("The following options are not permitted for drop.k: 'model','datafunction','inits', 'monitor', 'data', 'n.chains', 'combine', 'datalist', 'initlist', 'export.cluster', 'test'")
	}
	
	# Harvest required info:
	model <- runjags.object$model
	n.chains <- length(runjags.object$end.state)
	
	# These will be used quite a lot - arrays are explicitly rolled out to make overwriting easier:
	basedata <- runjags:::getjagsnames(list.format(runjags.object$data))
	testinits <- lapply(runjags.object$end.state, function(x) return(runjags:::getjagsnames(list.format(x))))
	baseinits <- lapply(runjags.object$end.state, list.format)
	
	getarraynames <- runjags:::getarraynames
	getjagsnames <- runjags:::getjagsnames
	
	# Match the data to vars:
	jagsnames <- names(basedata)
	vars <- jagsnames[runjags:::matchvars(runjags:::checkvalidmonitorname(dropvars), jagsnames)]
	newinits <- basedata[vars]
	# Remove the baseinits that overlap with the dropvars variable (data that is already missing):
	fake <- numeric(length(vars))
	names(fake) <- vars
	torem <- names(runjags:::getarraynames(fake))
	torem <- torem[torem %in% names(baseinits[[1]])]
	if(length(torem)>0){
		## Give a warning here??
		baseinits <- lapply(baseinits, function(x) return(x[!names(x)%in%torem]))
	}
	
	if(k > length(vars)){
		stop('Specified "k" was greater than the number of data points', call.=FALSE)
	}
	if(k < 1)
		stop('Specified "k" was less than 1', call.=FALSE)

	# One way to specify drop all:
	if(k == length(vars)){
		if(!is.na(simulations))
			warning('Value for "simulations" ignored as k == length(vars)', call.=FALSE)	
		simulations <- 1
	}
	
	if(length(simulations)!=1 || is.na(simulations)){
		if(k!=1)
			stop('A single positive integer must be provided for "simulations" if k>=2', call.=FALSE)		
		simulations <- length(vars)
	}
	
	# Other way to specify drop all:
	if(simulations==1){
		if(k!=1 && k!=length(vars)){
			warning('Value for "k" ignored as "simulations" was specified as 1', call.=FALSE)
		}
		k <- length(vars)
	}
	
	# The inits function takes base inits and base data and data removed and acts appropriately
	
	# datafunction should be able to see basedata and baseinits and vars as this is the parent frame
	dfdone <- FALSE
	avs <- FALSE
	if(k==1 && simulations==length(vars)){
		# Non random drop-1:
		datafunction <- function(i){
			thisdata <- basedata
			thisdata[vars[i]] <- as.numeric(NA)
			return(runjags:::getarraynames(thisdata))
		}
		dfdone <- TRUE
		avs <- FALSE
	}
	
  if(k==length(vars)){
		# Non random drop-all - this is probably not a useful thing to want to do....
		datafunction <- function(i){
			thisdata <- basedata
			thisdata[] <- as.numeric(NA)
			return(runjags:::getarraynames(thisdata))
		}
		dfdone <- TRUE
		avs <- FALSE
		if(simulations!=1 && !is.na(simulations))
			warning('Ignoring provided value for simulations, as "k" is equal to the number of variables')
		simulations <- 1		
	}
	if(!dfdone){
		# Otherwise, random drop-k:
		datafunction <- function(i){
			thisdata <- basedata
			makena <- sample(1:length(vars), k)
			thisdata[vars[makena]] <- as.numeric(NA)
			return(runjags:::getarraynames(thisdata))
		}		
		avs <- TRUE
	}
	
    inits <- function(chain){
		# We have the data returned by datafunction visible:
		thisdata <- runjags:::getjagsnames(data)
		# And the baseinits exported onto the cluster:
		baseinits <- baseinits
		# And the basedata exported onto the cluster:
    	basedata <- basedata
		# And the new inits:
		# Taking them from newinits only works if the array is complete - safer to take everything in data and make it NA:
		# thisinits <- newinits
		# makena <- which(!is.na(thisdata[names(thisinits)]))
		# thisinits[makena] <- as.numeric(NA)
		thisinits <- thisdata
		thisinits[] <- as.numeric(NA)
		# Find which are NA in the data:
		bringback <- names(thisdata)[is.na(thisdata)]
		thisinits[bringback] <- basedata[bringback]
		# And return all inits:
		toreturn <- c(baseinits[[chain]], runjags:::getarraynames(thisinits))
		toreturn <- toreturn[sapply(toreturn, function(x) return(!all(is.na(x))))]
		# Needs to be dump format to suppress finding #inits# from the model on the cluster:
		return(dump.format(toreturn))
    }

	results <- run.jags.study(simulations, model=model, datafunction=datafunction, data='', targets=basedata[vars], n.chains=n.chains, inits=inits, export.cluster=list(baseinits=baseinits, basedata=basedata, newinits=newinits, getjagsnames=getjagsnames, getarraynames=getarraynames), forceaverage=avs, ...)
	
	return(results)
	
}
