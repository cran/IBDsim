.findCAFs <-
function(x, ones, twos, loops = NULL) {  #common ancestral founders
    caf <- fou <- x$founders 
	fou_one = intersect(fou, ones)
	if (length(fou_one) > 1 || any(twos %in% fou)) return(numeric()) #since founders always have different alleles (by definition)
	
	fou_desc = lapply(fou, descendants, x=x, original.id=FALSE)
	
	if(length(twos)>0) {
		.ancestral.founders = function(id)	if (id %in% fou) return(id) else fou[sapply(fou_desc, function(ds) id %in% ds)]
		if(is.null(loops)) loops = pedigreeLoops(x)
		if(length(loops) == 0) return(numeric())  #NB: 1)pedigreeLoops bruker originale ids 2) veldig ineffektivt (men selv telemark går noen sekunder)
		bottoms = sapply(loops, '[[', 'bottom') 
		tops = sapply(loops, '[[', 'top')
		for (bot in twos) caf = intersect(caf, unlist(lapply(tops[bottoms==bot], .ancestral.founders)))
	}
	if(length(fou_one)==1) {
		test_ones = all(ones[ones != fou_one] %in% fou_desc[[match(fou_one, fou)]])
		if(test_ones && fou_one %in% caf) 	return(fou_one) 
		else return(numeric())
	} 
	else	
		caf[ sapply(caf, function(f) all(ones %in% fou_desc[[match(f, fou)]])) ]  #keep only cafs that have all the 'ones' among their descendants.
}


obligate.carriers = function(x, sap) {
	one = c(sap[['1']], sap[['atleast1']]); two = sap[['2']]
	loops = pedigreeLoops(x)
	bottoms = sapply(loops, '[[', 'bottom')
	tops = sapply(loops, '[[', 'top')
	if(any(one %in% bottoms)) stop("Inbred individual carrying 1 allele: Not implemented yet.")
	singleloop=function(inb, caf) {
		loop = loops[(bottoms == inb) & (tops == caf)]
		if(length(loop) != 1 ) stop("No loop, or more than one loop. Problem!")	
		loop=loop[[1]]
		c(loop[['top']], loop[['pathA']], loop[['pathB']])
	}
	cafs = .findCAFs(x, one, two, loops)
	caf_descpaths = paramlink:::.descentPaths(x, cafs, original.id=FALSE)
	
	#---non-inbred 'one' individuals: For each CAF there are no choices here---
	noninb_ones = setdiff(one, bottoms)
	obligate = lapply(cafs, function(caf) {
		paths = caf_descpaths[[match(caf, cafs)]]
		obl = sap[['1']] #initialize.
		for(id in noninb_ones) obl = c(obl, unlist(lapply(paths, function(pth) pth[seq_len(match(id, pth, nomatch=1)-1)])))
		for (inb in two) obl = c(obl, singleloop(inb, caf))
		sort.default(unique(obl))
	})
	setdiff(obligate, sap[['atleast1']])  #return list of vectors: these are the sample space for the '1'-element of the condition sap.
}
	

.whichIndivsIBD <-
function(x, share) {
	if (is.null(x$model)) chrom="AUTOSOMAL" else chrom=x$model$chrom
	if (is.null(x$sim)) sim=rep.int(2, x$nInd) else sim=x$sim
	
	if (is.character(share) && share %in% c("sim", "disease"))
		switch(chrom, 
		 "AUTOSOMAL" = {share = which(sim==2 & x$pedigree[,'AFF']==2); nonshare = which(sim==2 & x$pedigree[,'AFF']==1)},
		 "X" = {share = which(sim==2 & x$pedigree[,'SEX']==1 & x$pedigree[,'AFF']==2); nonshare = which(sim==2 & x$pedigree[,'SEX']==1 & x$pedigree[,'AFF']==1)}
		)
	else if (is.list(share) & length(share)==2) {share.list = share; share=share.list[[1]]; nonshare=share.list[[2]]}
	else if (is.numeric(share)) nonshare=numeric()
	return(list(yes=share, no=nonshare))
}

	
	
	
inbreeding = function(x) {
	ped = x$pedigree
	kin.coeff = kinship(id=ped[,'ID'], father.id=ped[,'FID'], mother.id=ped[,'MID'])
	inb.coeff = numeric()
	inb.coeff[x$founders] = 0
	inb.coeff[x$nonfounders] = sapply(x$nonfounders, function(i) kin.coeff[ped[i, 'FID'], ped[i, 'MID']])
	inb.coeff
}
