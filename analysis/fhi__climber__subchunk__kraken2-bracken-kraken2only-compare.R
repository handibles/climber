
rm(list=ls())

##  kraken2only  ===================================================================

# deleted # from top of file first!
krak2 <- read.table("output/r0937/krakenStnd_kreportCombo_count.tsv", header=TRUE, sep="\t", quote = "")
rownames(krak2) <- paste0("k2_", stringr::str_pad(krak2$taxid, width=7, pad=0))
head(krak2)
# View(krak2)

# choose to keep the abundances at "lvl", rather than at "tot", expect that 
kfeat <- krak2[ , grep("\\d_lvl$", colnames(krak2), value= TRUE, perl=TRUE) ]
head(kfeat)

# have to be smart and combine de novo, other tax doc doesnt have the S1 to help
tax_x <- data.frame(taxon=krak2[ ,30], k2_id=krak2[ ,29], row.names = rownames(krak2) , stringsAsFactors = FALSE)
tax_x[,1] <- gsub('\\/','', tax_x[,1])
tax_x[,1] <- gsub('ALOs_','ALOs-', tax_x[,1])
tax_x[,1] <- gsub('[\\[\\]]','', tax_x[,1], perl = TRUE)
tax_x$G <- gsub("^\\s*(\\w*)\\s.*", "\\1", tax_x$taxon, perl=TRUE)
tax_x$S <- gsub("^\\s*(\\w+)\\s(\\w+)\\s*.*", "\\2", tax_x$taxon, perl=TRUE)
head(tax_x)
# View(tax_x)

# now combine by "Genus species"
uniq <- apply( unique(tax_x[,c("G", "S")]), 1, function(aa){ rownames(tax_x)[(tax_x$G == unlist(aa[1])) & (tax_x$S == unlist(aa[2]))] })    #  aa <- unique(tax_x[,c("G", "S")])[240,]
kfeat_g <- t(sapply( uniq, function(aa){ colSums(kfeat[ aa , ]) }))

## < !! > not solved
poss_tax <- apply( unique(tax_x[,c("G", "S")]), 1, function(aa){ tax_c[ grepl( unlist(aa[1]), tax_c$G) & grepl( unlist(aa[2]), tax_c$S) , ] })  # aa <- unique(tax_x[,c("G", "S")])[230,]

##   standard kraken2+bracken  =====================================================================

krak2_brack_file="output/r0937/krakenStnd_abundances.tsv"
krak2_mpa_file="output/r0937/krakenStnd_taxonomy.tsv"

k2dat <- read.table(file = krak2_brack_file, header=TRUE, quote = "", sep = '\t')
k2dat[1:10, 1:10]

prenames_k2dat <- paste0("k2_", stringr::str_pad(k2dat[,2], width=7, pad=0))
## non-unique names -  assume only one dupe of each
prenames_k2dat[ duplicated(k2dat[,2])] <- gsub('k2_', 'k2-A_', prenames_k2dat[ duplicated(k2dat[,2])])
rownames(k2dat) <- prenames_k2dat

# counts
k2_counts <- k2dat[ , grep('_num', colnames(k2dat))]
colnames(k2_counts) <- gsub( ".bracken_num", "", colnames(k2_counts) )


##  taxa ===============================

tax_a <- data.frame(taxon=k2dat[ ,1], k2_id=k2dat[ ,2], row.names = rownames(k2dat) , stringsAsFactors = FALSE)
tax_a[,1] <- gsub('\\/','', tax_a[,1])
tax_a[,1] <- gsub('ALOs_','ALOs-', tax_a[,1])
head(tax_a)
tax_b <- read.table( krak2_mpa_file, sep='\t', header=FALSE, fill = TRUE, stringsAsFactors = FALSE, quote="")
head(tax_b)

## regularise taxonomic ranks (not all have K:P:C:O:F:G:S etc.)   -----------------
ranks <- c("k__", "p__", "c__", "o__" , "f__", "g__", "s__")
tax_b <- t(apply(tax_b, 1, function(aa){  # aa <- tax_b[1,]
  sapply(ranks, function(aaa){ # aaa <- "c__"
    bbb <- grep(aaa, unlist(aa), value=TRUE)
    ifelse( length(bbb) == 0, "unkn.", bbb)
  })
}))
dim(tax_b) ; str(tax_b)

## more tax issues due to names: doubled ranks:
prior_index <- (unlist( lapply(1:nrow(tax_b), function(aa){ if( sum( grepl("unkn.", tax_b[aa,])) == 6){aa} }) ) - 1)
tax_b[ prior_index , "s__"] <- tax_b[ (prior_index+1) , "s__"]
tax_b <- tax_b[ !apply(tax_b, 1, function(aa) sum(grepl("unkn.", unlist(aa)) ) == 6) , ]
colnames(tax_b) <- c("D", "P", "C", "O", "F", "G", "S")


tax_a[ , "taxon"] <- gsub(" ", "_", paste0("s__", tax_a[ , "taxon"]))
head(tax_a) ; head(tax_a[ , "taxon"]) ;  dim(tax_a)
head(tax_b) ;   head(tax_b[ , "S"]) ;  dim(tax_b)


## stitch  &  name   ----------------------------------

tax_c <- merge(tax_b, tax_a, by.x="S", by.y="taxon")[, c(2:7,1,8)]
head(tax_c) ; dim(tax_c)
prenames_tax <- paste0("k2_", stringr::str_pad(tax_c[,"k2_id"], width=7, pad=0))
rownames(tax_c) <- prenames_tax


## check all on the same page   -----------------------

all(rownames(k2dat) == rownames(k2_counts))
all(rownames(k2dat) %in% rownames(tax_c))

kraken_ids <- sort(rownames(k2dat))
k2dat <- k2dat [kraken_ids , ]
k2_counts <- k2_counts [kraken_ids , ]
tax_c <- tax_c[kraken_ids , ]

## remove if no values (kraken2 set to give all taxa)
pos_abund <- rownames(k2_counts)[rowSums(k2_counts) > 0]
k2_counts <- k2_counts[ pos_abund , ]
k2_tax <- tax_c[ pos_abund , ]


## check in vegan =================
library(vegan)

brack_g <- k2_counts[ , -1 ]
krak2_g <- kfeat[ rowSums(kfeat, na.rm=TRUE) > 10 , ]

# unify names and orders
colnames(krak2_g) <- gsub("X(\\d*)_.*", "nur_\\1", colnames(krak2_g), perl=TRUE)
colnames(brack_g) <- gsub("X_N(\\d*)_.*", "nur_\\1", colnames(brack_g), perl=TRUE)
brack_gs <- brack_g[ , sort(colnames(brack_g)) ]
krak2_gs <- krak2_g[ , sort(colnames(krak2_g)) ]

brack_gs[1:10,1:10]
krak2_gs[1:10,1:10]

str(brack_gs)
str(krak2_gs)

plot( colSums(brack_gs), colSums(krak2_gs)) 

brack_gRA <- t(apply(brack_g, 2, function(aa) aa/sum(aa)))
krak2_gRA <- t(apply(krak2_g, 2, function(aa) aa/sum(aa)))

brack_bc <- vegdist( brack_gRA, method = "bray")
krak2_bc <- vegan::vegdist( krak2_gRA, method = "bray")

par(mfrow = c(1,3))
#
plot( wcmdscale(brack_bc), type = "none", main = "shotgun using kraken2+bracken\nwith default db, 15:3:10 params")
text(wcmdscale(brack_bc), display = "sites",cex = 2, col = "red")
#
plot( wcmdscale(krak2_bc), type = "none", main = "shotgun using kraken2 only\nwith default db, 15:3:10 params")
text(wcmdscale(krak2_bc), display = "sites",cex = 2, col = "red")

# relate the two distance matrices
as.matrix(krak2_bc)[1:10, 1:10]
as.matrix(brack_bc)[1:10, 1:10]
(proc <- procrustes(brack_bc, krak2_bc))
# summary(proc)
plot(proc, main = "data are highly correlated (r=0.97, p=0.001)", kind=1)
(mant <- mantel(brack_bc, krak2_bc))


