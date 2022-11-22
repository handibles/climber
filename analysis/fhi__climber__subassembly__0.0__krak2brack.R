# move it somewhere sensible
mfile="input/fhi_redch_krakenStnd_abundances.tsv"

k2dat <- read.table(file = mfile, header=TRUE, quote = "", sep = '\t')
k2dat[1:10, 1:10]

prenames_k2dat <- paste0("k2_", stringr::str_pad(k2dat[,2], width=7, pad=0))
## non-unique names -  assume only one dupe of each
prenames_k2dat[ duplicated(k2dat[,2])] <- gsub('k2_', 'k2-A_', prenames_k2dat[ duplicated(k2dat[,2])])
rownames(k2dat) <- prenames_k2dat

# counts
k2_counts <- k2dat[ , grep('_num', colnames(k2dat))]
colnames(k2_counts) <- gsub( ".bracken_num", "", colnames(k2_counts) )

# relab%
k2_relab <- k2dat[ , grep('_frac', colnames(k2dat))]
colnames(k2_relab) <- gsub( ".bracken_frac", "", colnames(k2_relab) )

tax_a <- data.frame(taxon=k2dat[ ,1], k2_id=k2dat[ ,2], row.names = rownames(k2dat) , stringsAsFactors = FALSE)
tax_a[,1] <- gsub('\\/','', tax_a[,1])
tax_a[,1] <- gsub('ALOs_','ALOs-', tax_a[,1])
# tax_a[,1] <- gsub("'","", tax_a[,1])      # none?

## NOTE use quote="" to avoid nightmare of special characters like '"`/,| etc in names
tax_b <- read.table('input/fhi_redch_krakenStnd_taxonomy.tsv', sep='\t', header=FALSE, fill = TRUE, stringsAsFactors = FALSE, quote="")
head(tax_b)

## regularise taxonomic ranks (not all have K:P:C:O:F:G:S etc.)
ranks <- c("k__", "p__", "c__", "o__" , "f__", "g__", "s__")
tax_b <- t(apply(tax_b, 1, function(aa){  # aa <- tax_b[1,]Ó
  sapply(ranks, function(aaa){ # aaa <- "c__"
    bbb <- grep(aaa, unlist(aa), value=TRUE)
    ifelse( length(bbb) == 0, "unkn.", bbb)
  })
}))
dim(tax_b) ; str(tax_b)

## e.g. eukarya have loads of ranks, so many empty sapces to be filled, e.g.
# tax_b[4000:4200,]

## more tax issues due to names: doubled ranks:
prior_index <- (unlist( lapply(1:nrow(tax_b), function(aa){ if( sum( grepl("unkn.", tax_b[aa,])) == 6){aa} }) ) - 1)
# # compare:
# tax_b[ prior_index , ]
# tax_b[ prior_index+1 , ]
# amend, by doubling up (will remove duplicates later)
tax_b[ prior_index , "s__"] <- tax_b[ (prior_index+1) , "s__"]
# row is only useful if has a species: cannot have a species, and hit 6, And be worth keeping
tax_b <- tax_b[ !apply(tax_b, 1, function(aa) sum(grepl("unkn.", unlist(aa)) ) == 6) , ]

colnames(tax_b) <- c("D", "P", "C", "O", "F", "G", "S")


tax_a[ , "taxon"] <- gsub(" ", "_", paste0("s__", tax_a[ , "taxon"]))
head(tax_a) ; head(tax_a[ , "taxon"]) ;  dim(tax_a)
head(tax_b) ;   head(tax_b[ , "S"]) ;  dim(tax_b)


## stitch  &  name
tax_c <- merge(tax_b, tax_a, by.x="S", by.y="taxon")[, c(2:7,1,8)]
head(tax_c) ; dim(tax_c)
prenames_tax <- paste0("k2_", stringr::str_pad(tax_c[,"k2_id"], width=7, pad=0))
rownames(tax_c) <- prenames_tax


## check all on the same page
all(rownames(k2dat) == rownames(k2_counts) & rownames(k2dat) == rownames(k2_relab))
all(rownames(k2dat) %in% rownames(tax_c))

kraken_ids <- sort(rownames(k2dat))
k2dat <- k2dat [kraken_ids , ]
k2_counts <- k2_counts [kraken_ids , ]
k2_relab <- k2_relab [kraken_ids , ]
tax_c <- tax_c[kraken_ids , ]


## remove if no values (kraken2 set to give all taxa)
pos_abund <- rownames(k2_relab)[rowSums(k2_relab) > 0]
k2_relab <- k2_relab[ pos_abund , ]
k2_counts <- k2_counts[ pos_abund , ]
tax <- tax_c[ pos_abund , ]


mgdat <- data.frame("sample" = colnames(k2_counts),
                    "type" = gsub("(.).*", "\\1", colnames(k2_counts), perl = TRUE),
                    "seqdepth" = colSums(k2_counts) )

## check - do all match?
dim(k2_counts)
dim(k2_relab)
dim(tax)
dim(mgdat)

# look
k2_counts[1:10, 1:10]

# View(k2_counts)
# hist( log10(k2_counts ))

## tsv output
write.table(k2_counts, "output/fhi__redch__krakenStnd-20__num..tsv")
write.table(k2_relab, "output/fhi__redch__krakenStnd-20__frac..tsv")
write.table(tax, "output/fhi__redch__krakenStnd-20__tax..tsv")
write.table(mgdat, "output/fhi__redch__krakenStnd-20__dat..tsv")

## RDS for R output - note samples are ROWS
saveRDS( t(k2_counts), "output/fhi__redch__krakenStnd-20__num.RDS")
saveRDS( t(k2_relab), "output/fhi__redch__krakenStnd-20__frac.RDS")
saveRDS( tax, "output/fhi__redch__krakenStnd-20__tax.RDS")
saveRDS( mgdat, "output/fhi__redch__krakenStnd-20__dat.RDS")
