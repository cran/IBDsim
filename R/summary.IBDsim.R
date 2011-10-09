summary.ibd <-
function(runs, merged=TRUE) {
	total = attr(runs, 'total_map_length_Mb')
	disease = all(sapply(runs,function(r) any(r[,'disreg']==1)))

	res <- sapply(runs, function(run) {
		if(merged) r = mergeRuns(run) else r = run
		count = nrow(r)
		all.lengths = as.numeric(r[,'end']-r[,'start'])
		len = sum(all.lengths)
		aver = ifelse(count>0, len/count, 0)
		stats = c(count.all = count, fraction.all = len/total*100, average.all = aver)
		
		if(disease) {
			r = run[run[,'disreg']==0, , drop=F]
			if(merged) r = mergeRuns(r)
			count = nrow(r)
			rlengths = as.numeric(r[,'end']-r[,'start'])
			len = sum(rlengths)
			aver = ifelse(count>0, len/count, 0)
			stats = c(stats, count.rand = count, fraction.rand = len/total*100, average.rand = aver)
		
			dis = run[,'disreg']==1
			if(sum(dis)!=1) {print(run); warning("More or less than one disease region")}
			else {
				disrun = run[dis, ]
				dis.length = as.numeric(disrun['end']-disrun['start'])
				stats = c(stats, dis.length = dis.length, dis.rank =rank(-c(dis.length, all.lengths), ties="first")[1])
			}
		}
		stats
	})
	no.seg = mean(res['count.all',]==0)
	summary = round(c(rowMeans(res), 'zero.segments(%)' = no.seg*100), 3)
	print(summary)
	invisible(res)
}


	# percent.count = quantile(res['count.random',], c(0, .01, .05, .10, .25, .5, .75, .9, .95, 1))
	# percent.totlength = quantile(res['total.fraction.random',], c(0, .01, .05, .10, .25, .5, .75, .9, .95, 1))
	
	# navn = c(" total_fraction(%)", "mean_length(Mb)", "count")
	# #cat("\nAll region averages:\n"); print(means[1:4]); 
	# cat("\nRandom region averages:\n");	vec = means[1:3]; names(vec)=navn; print(vec)
	# cat("\nMerged random region averages:\n"); vec = means[4:6]; names(vec)=navn; print(vec)
	# #cat("\nEstimated percentiles for length (Mb) of disease region:\n")
	# #print(round(quantile(disrun[, 'length'], c(0, .01, .05, .10, .25, .5, .75, .9, .95, 1)),2))
	# cat("\nAverage percentiles for length (Mb) of individual random IBD regions:\n")
	# print(round(rowMeans(percent.length),2))
	# cat("\nPercentiles for number of of random IBD regions per genome:\n")
	# print(percent.count)
	# cat("\nPercentiles for total fraction random IBD:\n")
	# print(round(percent.totlength,2))
	# res
	# runstats

	
#merges overlapping and adjacent segments
mergeRuns <- function(r) {
	if (nrow(r)<2) return(r)
	r = r[order(r[,'chrom'], r[,'start'], r[,'end']),]
	mergelist = lapply(unique(r[, 'chrom']), function(i) {
		m = rangeUnion(list(r[ r[, 'chrom']==i, 2:3]))
		m = cbind(i, m)
		colnames(m) = c('chrom', 'start', 'end')
		m
	})
	do.call(rbind, mergelist)
}
	

segment.length.percentile = function(runs, what="all")
		sapply(runs, function(r) r)
		
		# if (nrow(r)==0) return(rep.int(0,10)) 
		# quantile(r[,'length'], c(0, .01, .05, .10, .25, .5, .75, .9, .95, 1))
		#
