
## TOTALLY NOT FINISHED AND NOT WORKING

# need to incorporate outptu and rank control

# assumes that set control carried out via sample names - and thats horseshit!


venn_ranks <- function(phylo = phylo , set = set , output = c('plot' , 'ranks', 'both'), ranks = c(1:6)){
  
  require(phyloseq)
  
  # <!> assumes set attribute is in otu_colnames
  
  # get the set compositions    <!> WARN :: assumes set attribute is in otu_colnames    <!> 
  x <- lapply(levels(sample_data(wmi)$month) , function(r) otu_table(wmi)[ , grep( r, colnames(otu_table(wmi)))] )  
  x <- sapply(x , function(q){ rownames(q[ rowSums(q) > 0 , ]) })
  names(x) <- levels(sample_data(wmi)$month)
  
  # use gplots::venn to get plotting
  gplots::venn(x) 
  
  # use gplots::venn to get intersections
  z <- gplots::venn(x, show.plot = FALSE)
  z <- attributes(z)$intersections
  y <- lapply(z, function(p) tax_table(wmi)[p , 1:6])
  names(y) <- names(z)

  # write sets and all intersections to ../output
  dir_out <- paste0('../output/venn_intersects_','wmi','\\.','month')
  dir.create(dir_out)
  lapply(names(y), function(w) write.table(y[w], paste0(dir_out , '/venn_overlap_set_', w,'_table.txt'), sep='\t'))

  # charts
  pdf(paste0(dir_out , '/venn_intersections_', 'month', '_plots.pdf'))
  lapply(names(y), function(w){
    a <- prune_taxa( (taxa_names(wmi_5RA) %in% rownames(y[w][[1]])), wmi_5RA)
      # a <- subset_samples(a, month %in% b )
        # b <- unlist(strsplit(w, split = ':'))
        # sample_data(wmi)$month %in% b }) # not crucial as should be subset out, but infuriating
    plot_bar(a, x='ID', fill='Genus') +
      facet_grid(Phylum~month, space = 'free', scales = 'free_x') +
      scale_colour_manual(values = e.col) +
      theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
      ggtitle(w)
    })
  dev.off()
  
  
  return(y)
}

venn_ranks(phylo = wmi, set = part, output = 'both')
