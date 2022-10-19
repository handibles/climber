## phylo_out

## README  ::  use like:
# phylo_out('~/Downloads/merge_phyloseq.RDS' , output_name = 'merge_ITS' , id='FYP')


phylo_out <- function(phylo , output_name = 'phylo_out' , out=NULL, id='ASV'){

  # ## testing
  # phylo <- '~/Downloads/merge_phyloseq.RDS'
  # id <- 'FYP_ASV'
  # output_name <- 'fyp_phylo'
  # out <- NULL

  p <- NULL
  if (class(phylo) == 'phyloseq'){ p <- phylo } else { p <- readRDS(phylo)}
  if (is.null(p)){ print('couldnt dig your phylo vibe my dogge; try submitting a *valid* phyloseq-class object, or the *valid* location of a phyloseq object saved as RDS')}
  
  library(stringr)

  # config - names, ID format etc
  pad_width <- length(unlist(strsplit(as.character(ntaxa(p)) , split = NULL)))
  uniqueseqs <- taxa_names(p)
  taxa_names(p) <- paste0(id, '_' , str_pad(1:ntaxa(p) , pad_width, pad = 0))
  
  # add names to tax_table, especially handy if original names are seqs
  t <- data.frame(tax_table(p), stringsAsFactors = FALSE)
  # rownames(t) <- taxa_names(p)
  t[,id] <- uniqueseqs
  tax_table(p) <- tax_table(as.matrix(t))

  if (is.null(out)){dir <- './'} else {dir <- paste0(out , '/')}
    
  write.table(t(otu_table(p)) , paste0(dir , output_name, '_otu_table.tsv') , sep='\t')
  write.table(tax_table(p) , paste0(dir , output_name, '_tax_table.tsv') , sep='\t')
  write.table(sample_data(p) , paste0(dir , output_name, '_sample_data.tsv') , sep='\t')
  
  # # crappy FNA out, assumes ASV in row 8 ; incurs terrible printout via cat
  # sapply(taxa_names(p), function(x){
  #   cat(paste0('>' , x), tax_table(p)[ x , 8], '',  file = paste0(dir , output_name, '_seqs.fna'), sep = '\n',append = TRUE)
  # })
  
  }





## jfg_21.03.19
