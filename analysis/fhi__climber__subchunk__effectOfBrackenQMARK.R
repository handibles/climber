## does bracken make a difference?

# stnd 
krak2_brack_file="input/fhi__redch_krakenStnd.0.15_abundances.tsv"
##
k2dat <- read.table(file = krak2_brack_file, header=TRUE, quote = "", sep = '\t')
prenames_k2dat <- paste0("k2_", stringr::str_pad(k2dat[,2], width=7, pad=0))
## non-unique names -  assume only one dupe of each
prenames_k2dat[ duplicated(k2dat[,2])] <- gsub('k2_', 'k2-A_', prenames_k2dat[ duplicated(k2dat[,2])])
rownames(k2dat) <- prenames_k2dat
# counts
k2_counts <- k2dat[ , grep('_num', colnames(k2dat))]
colnames(k2_counts) <- gsub( ".bracken_num", "", colnames(k2_counts) )
k2_counts[1:10,1:10]


## sh for kraken
grep -E \ts\t $KROUT/${TEST}_test_kraken_report ////