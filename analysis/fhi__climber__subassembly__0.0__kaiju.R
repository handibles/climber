
mfile="input/fhi__redch__kaiju_speciesab.tsv"
# mfile="input/nur_tot_kai_spec.tsv"


## import data into R.    ------------------

# note the quote option, to stop strange names doing strange things
kaidat <- read.table(file = mfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
# remove 1+ header(s)
kaidat <- kaidat[ !grepl("file", kaidat[,1]) , ]
dim(kaidat)
# View(kaidat)

## fix sample names
kaidat$file <- gsub(".*/_(N.*)_Nuria.*", "\\1", kaidat$file)
# kaidat$file <- gsub(".*/(.*)_kaiju.*", "\\1", kaidat$file)
kaidat$otu <- gsub( ".*;(.*);.*;", "\\1", kaidat$taxon_name )
kaidat$ID <- paste0( "kai__", stringr::str_pad( kaidat[, "taxon_id"], width = 7, pad = 0) )
length(unique(kaidat$file))         # 20 files
length(unique(kaidat$otu))          # 1.3K OTUs
kaidat$percent <- as.numeric(as.character(kaidat$percent))
sum(kaidat$percent, na.rm = TRUE)

head(kaidat) # View(kaidat)
str(kaidat)


## assess total output    -----------------

kaidf <- as.data.frame(kaidat) 
kaidf[ , "reads"] <- as.numeric(kaidf[ , "reads"])

sum( (kaidf[ , "reads"]))
sum( (kaidf[ , "percent"]), na.rm = TRUE)  # nsample*100%

## nreads unclassified / sent to bin
kaidf$substrate <- gsub(".*-([A-Z]*).*_.*", "\\1", kaidf$file)
# ggplot(kaidf, aes(log10(reads), fill = substrate)) + geom_boxplot() + labs(title = "total N reads per substrate")


## congeal  -------------------------------------------

u_bio <- unique( kaidat$otu)
u_sample <- unique( kaidat$file)

# to each unique taxon found (u_bio), assign the number of counts found in each unique sample (u_sample)
kaigeal <- do.call("rbind", lapply( u_bio, function(aa){   # aa <- u_bio[100]
  bb_df <- kaidat[ grep(aa, kaidat$otu) , ]
  cc_vec <- unlist(lapply( u_sample, function(aaa){ # aaa <- u_sample[9]
    as.numeric(bb_df[ match(aaa, bb_df$file) , "reads" ])
  }))
} )
)

rownames(kaigeal) <- u_bio# [1:100]
colnames(kaigeal) <- u_sample
# only keep real numbers
kaigeal[ is.na(kaigeal)] <- 0
dim( kaigeal <- kaigeal[ , !apply( kaigeal, 2, function(aa){ all(aa == 0)}) ] )

# look
kaigeal[1:10,1:10]

## taxonomy   -----------------------------------------

tax <- unique( t(# do.call("rbind", 
  apply( kaidat[ , c("ID", "taxon_name")], 1, function(aa){   # aa <- kaidat[ 1, c("ID", "taxon_name" )]
    bb_vec <- aa[2] # paste0( aa, collapse = ";")
    if( grepl("Viruses", bb_vec)){ 
      cc_vec <- c( rep( "Viruses", length(2:7)), aa[1])
    }else{
      cc_vec <- c( unlist(strsplit(unlist(bb_vec), ";"))[ 2:7], aa[1] )
    }
    cc_vec  }) ))
rownames(tax) <- tax[,5]
colnames(tax) <- c("p", "c", "o", "f", "g", "s", "ID")   # Jun'22 - remove "d" here, increase index of all names by 1
head(tax)
dim(tax <- (unique(tax)))


