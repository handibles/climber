#!/usr/bin/env Rscript


## R E A D   P R O V E N A N C E


library('dplyr')
library('dada2')
library('stringr')
library("ShortRead")
library("Biostrings")
library('seqTools')


  ## bash bit 
  
  ## Q U E R Y   L E N G T H S
  #improve with awk probably (2 cols, sort) - needs to be robust          <   !   >
  
  ## file, total reads, read length*count  ::  see https://www.biostars.org/p/72433/ for awk trick
  for i in $(ls /data/jamie/eoe/Materials/0_raw_reads/*.fastq | cut -f 7 -d "/");
  do echo $i ;                       # name
     cat $i | grep -c '^+$' ;        # read count ::   ^ , $ = BOL/EOL
     cat $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c ; # unique read lengths, per 4-line piece
  done > ~/Materials/0_read_profile.txt


## IBD sample space - i.e. not subset into different scripts

## Functionable -  really important, incorporate.

# fn for unique sequences
getN <- function(x) sum(getUniques(x))    #from shortReads

# UniqSeq counts to df
  track <- cbind( sapply(cutFs, getN),  # replace 'out' matrix made through getN() on filterAndTrim() objects 
                  sapply(cutRs, getN),
                  sapply(filtFs, getN),
                  sapply(filtRs, getN),
                  sapply(dadaFs, getN),  
                  sapply(dadaRs, getN),
                  sapply(mergers, getN),
                  rowSums(seqtab.nochim)
  )
  track <- as.data.frame(track)  #motherfucker
  colnames(track) <- c("inputF", "inputR", "filteredF", "filteredR", "denoisedF", "denoisedR", "merged", "ASV")
  rownames(track) <- sample.names
  track_bup <- track
  head(track)
  write.table(track,'ibd48_uniques_provenance.txt', sep='\t')



## =====================================

## T R A C K
# of interest: samples, groups, read count 

# have already sample/state mapping:
states <- read.table('IBD_total.map.txt', header=TRUE, sep='\t', row.names = 1)
rownames(states) <- gsub('\\.' , '_', rownames(states))   # match names
track <- track_bup

(rownames(track) %in% rownames(states))
track$state <- (states[rownames(track),])
track$sample <- rownames(track)

z.rm <- melt(track)
# View(z.rm)
# establish max read counts for plotting:
max(track[,1:6])

ggplot(z.rm, aes(y=value, x=variable,  fill=as.factor(state) ,  colour=as.factor(state) )) +
  geom_boxplot(outlier.shape = NA, colour='black') + 
  geom_jitter(width=0.2, shape=21,  size=2.5, colour='black', aes(fill=as.factor(state) )) + #note factoring-by-source requires use of the aes arg
  scale_fill_manual( "state",values=c("#1e8bc3",  "#C30000",  "#ffb90f") ) +         # from RY flow 
  scale_color_manual("state",values=c("#1e8bc3",  "#C30000",  "#ffb90f") ) +
  scale_y_continuous(limits = c(0, 150000)) +
  labs(title='Read Counts Through DADA2 pipeline' ) + 
  xlab('DADA2 Step') + 
  ylab('number of reads')

## ===========================================

## from phyloseq kick-off


## work this out, or into a standard graphics pipeline

# ## S U M M A R Y   S T A T S  ##
# # phylum breakdown across whole study, per-reactor, and maybe even per sample...
# # not the same as the average % per sample, i.e. irrespective of library size.
# #
# ##   M A K E   M E   I N T O   A   F U N C T I O N
# ### take phylo object, take sample variable: cut phylo by variable, by rank
# ## A L S O   A   F U N C T I O N   F O R   O T U   I D s
# 
# eoe_phyglom <- tax_glom(eoe,taxrank='Phylum')  # one glom instance preserves order in later sub-sets
# 
# # sum_phy inherentl ordered by SV abundance
# sum_phy <- as.data.frame( round( (taxa_sums(eoe_phyglom)/sum(taxa_sums(eoe) ) *100), digits=4 ),  # grab rounded percentages
#                           row.names = as.character(tax_table(eoe_phyglom)[,2]) )                 # w. phyla as row names , cant set col.names?
# 
# # can illustrate taxa tables are intact based on sort(taxa_sums(z.isa/b/bes/ces)) etc.
# z.isa <- subset_samples(eoe_phyglom,Description=='1.ISA' )
# sum_phy[,2] <- round( (taxa_sums(z.isa)/sum(taxa_sums(z.isa) ) *100), digits=4 )  #keep NArm to keep comparable order between reactors
# z.isb <- subset_samples(eoe_phyglom,Description=='2.ISB' )
# sum_phy[,3] <- round( (taxa_sums(z.isb)/sum(taxa_sums(z.isb) ) *100), digits=4 )
# z.bes <- subset_samples(eoe_phyglom,Description=='3.BES' )
# sum_phy[,4] <- round( (taxa_sums(z.bes)/sum(taxa_sums(z.bes) ) *100), digits=4 )
# z.ces <- subset_samples(eoe_phyglom,Description=='4.CES' )
# sum_phy[,5] <- round( (taxa_sums(z.ces)/sum(taxa_sums(z.ces) ) *100), digits=4 )
# 
# colnames(sum_phy) <- c('TOT%','ISA%','ISB%','BES%','CES%') # cant seem to sort (order) properly either
# 
# ## add number of genera per phylum (overall) #OR just the 3-5 phyla you'll actually talk about..
# #length(get_taxa_unique(subset_taxa(eoe,Phylum=='Firmicutes'),'Genus')) # this will overlook NAs etc
# dim(tax_table(subset_taxa(eoe,Phylum=='Firmicutes')))
# dim(tax_table(subset_taxa(eoe,Phylum=='Euryarchaeota')))
# dim(tax_table(subset_taxa(eoe,Phylum=='Bacteroidetes')))
# 
# View(sum_phy)
# 
# ## Stats per eoe_01 / eoe_un01
# # total ASV_orig, total sequences_orig
# eoe_orig ; sum(sample_sums(eoe_orig)) ; mean(sample_sums(eoe_orig))
# # total ASV, total sequences
# eoe ; sum(sample_sums(eoe)) ; mean(sample_sums(eoe))
# # reads / sample
# sum_phy 
# # total un_ASV, total un_sequences
# eoe_01 ; sum(sample_sums(otu_table(eoe)[,rownames(tax_table(eoe_01))])) ; mean(sample_sums(otu_table(eoe)[,rownames(tax_table(eoe_01))]))
# # total un_ASV, total un_sequences
# eoe_un01 ; sum(sample_sums(otu_table(eoe)[,rownames(tax_table(eoe_un01))])) ; mean(sample_sums(otu_table(eoe)[,rownames(tax_table(eoe_un01))]))
# # Firm / Bact / Eury % range
# sum_phy # see above
# # #firm / #bact / #eury
# dim(tax_table(subset_taxa(eoe,Phylum=='Firmicutes')))
# dim(tax_table(subset_taxa(eoe,Phylum=='Euryarchaeota')))
# dim(tax_table(subset_taxa(eoe,Phylum=='Bacteroidetes')))
# 
# # NON- Firm / Bact / Eury % range
# c(100, 100, 100, 100, 100) - colSums(sum_phy[1:3,])
# get_taxa_unique(eoe,taxonomic.rank = 'Phylum') ; length(get_taxa_unique(eoe,taxonomic.rank = 'Phylum'))
# #subtract this from 2300:
# (dim(tax_table(subset_taxa(eoe,Phylum=='Firmicutes'))) ) + ( dim(tax_table(subset_taxa(eoe,Phylum=='Euryarchaeota'))) ) + ( dim(tax_table(subset_taxa(eoe,Phylum=='Bacteroidetes'))) )
# 2300 - 1662
# # 

