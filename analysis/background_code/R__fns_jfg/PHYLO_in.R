
# goes in the master ledger
# could also be a function, but for all the messing with basic sampe data

library(phyloseq)

# get in files (assume taxa are cols, DADA2 default)
a <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
t <- tax_table(taxa)
summary( colnames(a) == rownames(t) )  # sanity, ASVs equal between objects

## amend
  # names of ASVs, retain sequence in tax_table, amend taxranks
  colnames(t) <- c("Domain","Phylum","Class","Order","Family","Genus", "Species", "ASV")[1:ncol(t)]

  # send seqs to tax_table
  asv_seq <- row.names(t)
  t <- cbind(t , asv_seq)
  summary(rownames(t) == t[,7])  # sanity 2
  
  # rename ASVs 
  id = paste0("ASV_", stringr::str_pad(1:dim(a)[2], nchar(dim(a)[2]), pad=0))  # from write_FNA
  colnames(a) <- id
  rownames(t) <- id

  
## basic sample data  
  s <- rownames(a)
  # Factor = sapply(s, function(x) strsplit(x, '-')[[1]][1])   # split name to data
  s<- data.frame(ID = s, 
                   A = sapply(s, function(x) strsplit(x, '-')[[1]][1]),
                   B = sapply(s, function(x) strsplit(x, '-')[[1]][3]),
                   C = sapply(s, function(x) strsplit(x, '-')[[1]][2]),
                   row.names = s,
                   stringsAsFactors = FALSE)

 # # ensure variable correctly ordered
 # s$A <- factor(s$A, levels = c('A1' , 'A2' , 'A3', 'A4', 'A5'))
 # s$B <- factor(s$part, levels = c('B1' , 'B2' , 'B3', 'B4'))  
 

## make phylo
  a <- otu_table(t(a))
  t <- tax_table(t)
  s <- sample_data(s)
  PHYLO <- merge_phyloseq(a, s, t)
  saveRDS(PHYLO, './input/PHYLO_phylo.RDS')

