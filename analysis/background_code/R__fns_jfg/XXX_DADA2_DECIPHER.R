#!/usr/bin/env Rscript

### HONESTLY, BELONGS IN 7.3
# Method should be robust to conglomeration of multiple different 'runs' as provided for by the run_sample_list.txt etc
# dangerous to go alone: take one of these
  # [training](http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r132_March2018.RData)
  # [reference](http://www2.decipher.codes/SILVA-117_LSU_rRNA_Reference_Database.sqlite.zip)

#   #===============================================================================================================#



suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(phyloseq))   # make and remerge the phlyo table with otu & tax
suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(stringr))

# pre set
def_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

"Usage: 
do some taxa assignments dot com.R [options]
Description:   Use the new IDTAXA method to assign taxoonomy, with ML to avoid over-classification of new IDs
Options:
--inputRDS=<inputRDS>      [default: . ]    Path to RDS with merged seqtab
--outdir=<outdir>          [default: . ]    Path to output the tax file
--taxref=<taxref>          [default: . ]    Path to reference database
--strand=<strand>          [default: top]   top, bottom, or both strands (for mixed orientation reads)
--cpu=<cpu>                [default: NULL]  how many beasties to use: NULL (default) = use all available
--verbose=<verbose>        [default: TRUE]  verbose output, T/F
--ranks=<ranks>            [default: def_ranks] char vector of taxonomic levels i.e. desired colnames(tax_table)
" -> doc

opts <- docopt(doc)

# check opts
inputRDS       <- opts$inputRDS
outdir         <- opts$outdir
taxref         <- opts$taxref
strand         <- opts$strand
cpu            <- opts$cpu
ranks          <- opts$ranks

  # ## perennial jfg hack
  # inputRDS       <- '/claesson/jamie/bl_jamiedata/eoe/Materials/6_eoe_Seqtab/eoe_run1_seqtab.RDS'
  # outdir         <- '/claesson/jamie/bl_jamiedata/eoe/Materials/6_eoe_Seqtab/'
  # taxref         <- '/claesson/jamie/ref/SILVA_SSU_r132_March2018.RData'
  # strand         <- 'top'
  # cpu            <- 18
  # ranks          <- def_ranks
  # verbose        <- TRUE


# read in data from docopt.ions
seqtab.nochim <- readRDS(inputRDS)   # inelegant, but very transparent
trainingSet <- load(taxref)          # this probably incorrect


  # from DADA2 1.8
  dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
  ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose= verbose ) # use all processors
  ranks <- def_ranks # ranks of interest, set above docopts options

  # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
  taxid <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))

  # save names, give names
  taxid <- cbind(taxid,as.matrix(getSequences(seqtab.nochim)))
  colnames(taxid) <- c(ranks, "ASV")
  rownames(taxid) <- paste0("Seq_", str_pad(seq(nrow(taxid)), 7, pad = 0))
  
  rownames(tax_phy) == colnames(otu_table(phylo))
  
  # save output taxa as tsv and in phylobject
  tax_phy <- tax_table(taxid) ; saveRDS(tax_phy, 'eoe_taxtable.RDS')
  phylo <- readRDS('/claesson/jamie/bl_jamiedata/eoe/Materials/6_eoe_Seqtab/eoe_run1_phyloseq.RDS')
  tax_table(phylo) <- tax_phy    # using merge_/phyloseq doesn't work?
  phylo ; saveRDS(phylo , '/claesson/jamie/bl_jamiedata/eoe/Materials/6_eoe_Seqtab/eoe_run1_phyloseq.RDS')
  
  combineRDSfiles(RDSdata,
                  otutableout =     paste0( outdir,"/", experimentname, "_otus.txt"),
                  phyloseqoutfile = paste0( outdir,"/", experimentname, "_phyloseq.RDS"),
                  FNAoutfile =  paste0( outdir,"/", experimentname, "_otus.fna")
  )
  
  
q('no')