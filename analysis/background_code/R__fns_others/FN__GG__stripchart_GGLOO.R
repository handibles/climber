# a function to plot strip charts at various taxonomic levels, cutoffs and x axes
aldex.stripchart <- function(
  aldex.out=x.all, group.label="genus", 
  x.axis="diff.btw", cex=0.8, cutoff=0.1)
{
  # aldex.out holds the new data frame
  # group.label is the column name from the taxon file
  # x.axis should be either diff.btw (difference) or effect
  # cex controls plot size
  # cutoff is the BH corrected wilcoxon rank test statistic
  # copy the aldex output data frame
  aldex.out <- data.frame(aldex.out, taxon[rownames(x.all), group.label])
  colnames(aldex.out)[ncol(aldex.out)] <- group.label
  # some R housekeeping
  aldex.out <- droplevels(aldex.out)
  
  # get a vector of significant and non-significant rows
  cutoff <- 0.1 # change this to whatever you want
  non.sig <- aldex.out$wi.eBH >= cutoff
  sig.pos <- aldex.out$wi.eBH < cutoff & aldex.out$effect > 0
  sig.neg <- aldex.out$wi.eBH < cutoff & aldex.out$effect < 0
  
  # make a list of unique taxa in the dataset
  # there may be empty groups, ignore these for now
  groups <- unique(aldex.out$group.label)
  
  # generate a y axis plotting limit a bit larger than needed
  ylim<-c(length(groups) - (length(groups)+0.5), length(groups) + 0.5)
  
  x.axis <- x.axis # we find the effect or diff.btw columns to be most useful
  xlim = c(min(-1 * max(abs(aldex.out[,x.axis]))), max(aldex.out[,x.axis]))
  
  # basically can call the different significant groups within strip chart
  par(mar=c(5,15,5,1), las=1, cex=cex)
  stripchart(aldex.out[non.sig,x.axis] ~ aldex.out[non.sig,group.label],
             col=rgb(0,0,0,0.3), method="jitter", pch=19, xlim=xlim, xlab=x.axis,
             main=group.label)
  stripchart(aldex.out[sig.pos,x.axis] ~ aldex.out[sig.pos,group.label], 
             col=rgb(0,0,1,0.5),method="jitter", pch=19, add=T)
  stripchart(aldex.out[sig.neg,x.axis] ~ aldex.out[sig.neg,group.label], 
             col=rgb(1,0,0,0.5),method="jitter", pch=19, add=T)
  
  abline(v=0, lty=2, col=rgb(0,0,0,0.2),lwd=2)
}