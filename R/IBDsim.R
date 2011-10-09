IBDsim <-
function(x, sims, query=NULL, condition=NULL, map="decode", chromosomes=1:22, model="chi", merged=TRUE, simdata=NULL, skip.recomb = NULL) {
	if(!all(x$orig.ids==1:x$nInd)) stop("Individual ID's must be 1, 2, 3, ... . Please relabel (see e.g. ?relabel).")
	starttime = proc.time()
		
	if (is.null(simdata)) {
		map = loadMap(map, chrom = chromosomes)
	
		if(!is.null(condition)) {
			if(length(map)==1) dischr = rep.int(attr(map[[1]], 'chromosome'), sims)
			else dischr = sample(sapply(map, attr, 'chromosome'), size=sims, replace=T, prob=sapply(map, attr, 'length_Mb'))
			oblig.saps = sample.obligates(x, condition, sims)
		}
		else dischr=rep.int(0, sims)
	
		simdata = lapply(1:sims, function(i) 
			lapply(map, function(m) {
				if(dischr[sims]==attr(m, 'chromosome')) cond=oblig.saps[[i]] else cond=NULL 
				genedrop(x, map=m, condition=cond, model=model, skip.recomb=skip.recomb)
			}))
		attr(simdata, 'total_map_length_Mb') = attr(map, "length_Mb")
		cat("Simulation finished. Time used:", (proc.time()-starttime)[['elapsed']], "seconds\n")
	}
	if(is.null(query)) return(invisible(simdata))
	
	if(!is.null(two <- query[['2']]))
		{cat("\nInbreeding coefficients:\n"); print(cbind(ID=two, f=inbreeding(x)[two]))}
		
	runs <- lapply(simdata, function(h) 	sap.segments(h, sap=query))
	attr(runs, 'total_map_length_Mb') = attr(simdata, 'total_map_length_Mb')
	cat("\nResults:\n")
	summary.ibd(runs, merged=merged)
	cat("\nTotal time used:", (proc.time()-starttime)[['elapsed']], "seconds.\n")
	invisible(list(simdata=simdata, segments=runs))
}

sample.obligates = function(x, condition, sims) {
	obligate_ones = obligate.carriers(x, condition)
	complete.saps = lapply(obligate_ones, function(ones) {sap = condition; sap[['1']] = ones; sap})
	if(length(complete.saps)==1) {cat("Conditioning on the following SAP:\n"); .printSAP(complete.saps[[1]])}
	else {
		cat("Sampling condition SAPs among the following:\n\n") 
		for(i in 1:length(complete.saps)) {cat("SAP ",i,":\n",sep=""); .printSAP(complete.saps[[i]])}
	}
	weight = sapply(obligate_ones, function(vec) .5^(length(vec)-1))
	oblig.samples = sample(complete.saps, size=sims, replace=TRUE, prob=weight)
}


#Lag snarvei: condition = "AD" (alle syke bærer én), query = "AD" (bruker 'sim' til å lage sap'en)
	# condition = list(carry0 = ped[,'ID'][ped[,'AFF']==1], carry1 = ped[,'ID'][ped[,'AFF']==2])
	# search_pattern = list('0' = ped[,'ID'][x$sim == 2 & ped[,'AFF']==1], '1' = ped[,'ID'][x$sim == 2 & ped[,'AFF']==2]) 
	# #print(condition); print(search_pattern)
	

.printSAP = function(sap) {
	if(!is.null(two <- sap[['2']])) cat("Two copies:", paste(two, collapse=", "), "\n")
	if(!is.null(atl1 <- sap[['atleast1']])) cat("At least one copy:", paste(atl1, collapse=", "), "\n")
	if(!is.null(one <- sap[['1']])) cat("One copy:", paste(one, collapse=", "), "\n")
	if(!is.null(atm1 <- sap[['atmost1']])) cat("At most one copy:", paste(atm1, collapse=", "), "\n")
	if(!is.null(zero <- sap[['0']])) cat("Zero copies:", paste(zero, collapse=", "), "\n")
	cat("\n")
}

.prettycat = function(v, andor) {
	len <- length(v)
	switch(min(len, 3), toString(v), paste(v, collapse=" and "), paste(paste(v[-len], collapse=", "), andor, v[len]))
}