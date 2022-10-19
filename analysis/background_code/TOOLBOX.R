
# # # # # # # # # # # # # # # # # # # # #
# # # #     T O O L      B O X    # # # # 
# # # # # # # # # # # # # # # # # # # # #

# Check out GGally...


## jfg 4.11 tricks



## ========================================================================
    ##   W O R K I N G   L I N E T Y P E    P L O T
    # may it be a light for you in dark places

# ----------
eoe <- readRDS('./output/7_phyout/eoe_run1_phyloseq_uc_purMANUAL_GMPRnorm.RDS')
mgfeat <- as.matrix(data.frame(otu_table(eoe))) ; dim(mgfeat) 
mgdat <- data.frame(sample_data(eoe), stringsAsFactors = FALSE)
# ----------

mgfeat_ra <- t(apply(mgfeat, 1, function(x){ x/sum(x) }))

bc_dist <- vegan::vegdist(mgfeat_ra, method = "bray")
bc_pcoa <- ape::pcoa(bc_dist)
bc_df = data.frame(data.frame(mgdat, stringsAsFactors = FALSE),
                   "PCo1" = bc_pcoa$vectors[,1]*-1,
                   "PCo2" = bc_pcoa$vectors[,2]*-1)

ggplot(bc_df, aes(x = PCo1, y = PCo2)) + 
  stat_ellipse(aes(x = PCo1,y =PCo2, colour=Description, lty=Description),
               geom="polygon", alpha=0.0, show.legend = TRUE)
## =======================================================================



##   P R I M E  -   from howto parallel
library('parallel')
no_cores <- detectCores() - 3
cl <- makeCluster(no_cores)
# stopCluster(cl)



## for a feat_table j, metadata, and taxo, melt into columns

  ## bulk out last column with var of interest
  jjj <- melt(cbind(j_crush, j[, jdh$order]), variable_name = "Sample")    
  b <-  lapply( c("sample"), function(m){                # where sample is a feature name
    unlist( lapply( metadata[jdh$order, m], function(n) rep(n, nrow(taxo)) ))
  })

  ## or ACTUALLY melt & merge 
  #                            tax     feats
  jt <- melt(cbind.data.frame(  j_crush, j, stringsAsFactors=FALSE ))
  #               vars
  jtm <- (merge(jt, m, by.x="variable", by.y="sample"))


  # strip factors
  df <- as.data.frame(lapply(df, function(f){ if(class(f)=="factor"){as.character(f)}else{f}}), stringsAsFactors = FALSE)

  # comaprison of vlaues
    # axial plots
    #- display significance
    library('ggsignif')
    	+ geom_signif(comparisons = list(c("CONTROL", "POST", "PRE")), map_signif_level=TRUE)
  



## LOVELY subset, so simple
keep <- apply(counts, 2, function(x) sum(x >= 0.05) >= 10)


# list-of-lists to dataframe, beautifully, simply:
z <- recursive-list_structure
data.table::rbindlist(z)     #  THIS A THOUSAND TIMES THIS, think of the taxonomic possibilities

# with simple named lists, use base stack/unstack!
stack(z)


# arranging non-ggplots
par(mfrow=c(2,2))
plot(a)
plot(b)
plot(c)
plot(d)

z$i
# > Error in z$i : $ operator is invalid for atomic vectors  -  use *attributes* to access ye olde atomic vectors
attributes(z)$i
# > aww yeah

## upgrade bioconductor (including BiocManager?...)
BiocManager::install(version="3.X")

# reorder graph wihtout messing around with DFs and levels: do it wiht the actual es() call!
# https://stackoverflow.com/a/48998414

  level_order <- c('virginica', 'versicolor', 'setosa') #this vector might be useful for other plots/analyses
  ggplot(iris, aes(x = factor(Species, level = level_order), y = Petal.Width)) + geom_col()
  # or 
  level_order <- factor(iris$Species, level = c('virginica', 'versicolor', 'setosa'))
  ggplot(iris, aes(x = level_order, y = Petal.Width)) + geom_col()
  # or
  ggplot(iris, aes(x = factor(Species, level = c('virginica', 'versicolor', 'setosa')), y = Petal.Width)) + geom_col()

  
## split sample ID into meaningful data using base::srsplit()
# for df a, assuming delimited by '-'
s <- rownames(a)
Month = sapply(s, function(x) strsplit(x, '-')[[1]][1])


## RMD date
date: "`r format(Sys.time(), '%d %B, %Y')`"


## call dwn the moon
a <- c('dada2','phyloseq','DESeq2','ggplot2','seqTools','parallel','zCompositions','ALDEx2','compositions','knitr','igraph','phangorn','DECIPHER','ComplexHeatmap' ,'cowplot','shiny','shinythemes','tibble','networkD3','circlize','FSA')
sapply(a , function(x) install(x, ask = FALSE, update = TRUE))

## ~resolve missing packages
update.packages()


## find index of object (from vector X in vector Y)
sapply(X, function(x) match(x , Y))

## pad out names
paste0('fungi',formatC( 1:ntaxa(its), width=3,format='d',flag='0'))

## thank you hadley
# # It's recommended to use a named vector
# cols <- c("8" = "red", "4" = "blue", "6" = "darkgreen", "10" = "orange")
# p + scale_colour_manual(values = cols)

## bad melting, because categories numeric: bin them! 
bc_pcoa_boxp$years <- cut(bc_pcoa_boxp$years, c(-1,1,11,21,31,41), labels = c('never', '1-10', '11-20', '21-30', '>30'))


# this works for alphabetical sorting, but is dependent on what you're working with. if getting NAs, check factor matches between target/source
# with(tax_df, order(phylum, class, order))

# flatten motherfucking matrices  -  christ thank you Anna
unlist(lapply(tax_table(eoe_full)[,8] , as.character))
# or, far more simply, use data.frame (not AS.data.frame)
data.frame(tax_table(whatev), stringsAsFactors = FALSE)   # note stringsAsFactors = FALSE

#change axes (e.g. x axes) text: 
plot(..., xaxt = 'n') # first nuke original labels
axis(1, at=1:71, labels = names(sort(sample_sums(eoe))) , cex=0.6 , srt=90, line=4 ) #pos=1, offset = -3,  #then design your own, see SO via JD Long


# sexy tables:  RMD can format data you igve it for you, but if you want to format thge output from e.g. a df, try:
kable(df)  # from knitr (woo!), or
pander(df) # from pander, if you can get it to work. 

# split and index names (could gsub the rest out instead!)
samples <- sapply(fwd_reads, function(x) strsplit( x, '\\.')[[1]][2], USE.NAMES = FALSE )


# this will do colours for you between 3 and the max in its vector (9) 
brewer.pal(1, 'YlGn')
# this will do any number, including 1. Will also work on any number of specified colours, so no *need* to use brewer.pal
show_col(colorRampPalette(brewer.pal(10, 'YlGn'))(1))


# for each taxon in rank_order, order that taxons abundances and append it to df0
## either of these loops should work, but for presence of NAs (over- or under-matched)
df0 <- data.frame()
for (i in rank_order) {
  df0 <- rbind( df0 , (z.rm3[z.rm3$Class == i , ])[ with( z.rm3[ z.rm3$Class == i ,] , order(Abundance))  , ] )
}

df0 <- data.frame()
for (i in rank_order) {
  otu_order <- rownames(tax_table(gas_10)[ tax_table(gas_10)[ , 'Class'] == i , ])
  otu_order <- names(sort(taxa_sums(prune_taxa(otu_order , gas_10)), decreasing = TRUE))
  for (j in otu_order) {
    df0 <- rbind(df0 , (z.rm3)[ (z.rm3$OTU == j) , ] )
  }  
} 


# length-filter RDS sample data frames
for (samp in 1:length(RDSdata)) {          
  RDSdata[[samp]] <- RDSdata[[samp]][ apply(RDSdata[[samp]][1], 1, function(x) nchar(x) >= thresh) , ]  ##  says all reads _ARE_ 420 long
}

## filter based on 'run#' in filename
if (!(is.null(run))) {                                  # conditional for getting specific run
  RDSfiles <- grepl(paste0('run',run), RDSfiles)        ## C A R E F U L  , will nuke your shit if names are wrong     < ! >
  experimentname <- paste0(experimentname,'_run',run)   # append run# to output
}

## ugly but proud: make a list the length of number of read-pairs, in which to dump locations
readsin <- as.list(1:length(fnFs)) 
for(i in seq_along(fnFs)) {
  readsin[[i]] <- c( paste0(path, c(fnFs[i], fnRs[i]) ))  
}

## Not R, but hey  
## - X A R G S
# dumb copy
ls /data/sidney/16s/biopsy/run1/*gz | tail -n 500 | xargs cp -t /data/jamie/Materials/0_ry_reads
# dumb move
ls /data/jamie/Materials/3_d2 | tail -n 500 | xargs mv -t /data/jamie/Materials/3_d2/3_D2_REST

##  - T E E  :: Output to log
 "./aaa.sh 2>&1 | tee -a log" from
# https://stackoverflow.com/questions/692000/how-do-i-write-stderr-to-a-file-while-using-tee-with-a-pipe/692407#692407




# FRYAN/Sidney Scripts

# get variance per axes for PCA as follows:
z.pc <- prcomp(braycurtis_dist10)
biplot(z.pc)
 # Sum the total variance
 d.mvar <- sum(d.pcx$sdev^2)
 # get total `
 PC1 <- paste("GBM PC1: ", round(100*(sum(d.pcrx$sdev[1]^2)/d.mvar), 3), '%')   #making axes labels
 PC2 <- paste("GBM PC2: ", round(100*(sum(d.pcrx$sdev[2]^2)/d.mvar), 3), '%')


## for loop to find err_out's
fns = list.files()
dadaFs <- fns[grepl(".dadaF.RDS", fns)]
dadaF_list = list()
for (sample in dadaFs){
  x = readRDS(sample)
  dadaF_list[[sample]] = x$err_out       # insert err_out into list()'d object
}





## Sys and Env====== #

## install something from repo
source("http://bioconductor.org/biocLite.R")
devel = "http://bioconductor.org/packages/2.14/bioc"
biocLite("phyloseq", siteRepos=devel, suppressUpdates=TRUE, type="source")

## or from source:
install.packages(path_to_file, repos = NULL, type="source")
install.packages("https://bioconductor.org/packages/devel/bioc/src/contrib/ComplexHeatmap_1.99.5.tar.gz", repo=NULL, type="source")

## or from GH:
install_github("DeveloperName/PackageName")

# need to change Rscript library path, or keep re-setting location ::  https://stackoverflow.com/questions/27673000/rscript-there-is-no-package-called
.libPaths(c("/home/jfg/R/x86_64-pc-linux-gnu-library/3.5", .libPaths()))
# print off all the env variables
Sys.getenv()

# all the sourced libraries
.libPaths()
# [1] "/home/jfg/R/x86_64-pc-linux-gnu-library/3.5"
# [2] "/usr/local/lib/R/site-library"              
# [3] "/usr/lib/R/site-library"                    
# [4] "/usr/lib/R/library"                         

# ================== #


# # # # # # # # # # # # # # # # # # # # #
# # # #      P h D   T O O L S    # # # # 
# # # # # # # # # # # # # # # # # # # # #


## p h y l o    S U B S E T abund over 200 reads in at least 0.23 of samples
prunA = genefilter_sample(mv, filterfun_sample(function(x) x >=200), A=0.23*nsamples(mv)) 
mv_10 = prune_taxa(prunA, mv)
mv_un10 = prune_taxa(!prunA, mv)


# get numers in seq or sequence, or period, or step:
seq(0,200,20)

# export normaliZed counts from DESeq2
mv_in <- phyloseq_to_deseq2(mv, ~situ)
mv_in$situ <-relevel(mv_in$situ,ref='In-Situ')
mv_in_r<-DESeq(mv_in)
mv_D2 <- counts(mv_in_r,normalized=TRUE)

## phyloseq to DESeq
dw_d2 <- phyloseq_to_deseq2(dw_slv, ~Treat)
dw_d2 <- DESeq(dw_d2, test='Wald', fitType='parametric') 
## for Rlog or variance stabilising transformation (VST)
dw_d2_vst <- varianceStabilizingTransformation(dw_d2)
dw_d2_vst <- assay(dw_d2_vst)
write.table(dw_d2_vst,'dw_D2_vst.txt',sep='\t')

# avoid reverting to vector with single column dataframes: include drop=FALSE at end of index call
df[df$column == 'VALUE', , drop=FALSE]

# group plots together
grid.arrange(nrow = 1, plot1, plot2, plot3)

# ggplot's default palette, an even spacing of the hues of the colour wheel to a  required number (e.g. 10)
show_col(hue_pal()(10))
# good viridis palette
z.rm <- c('#21908CFF','#91908CFF','#27AD81FF') ; show_col(z.rm)

## combine two DF elements into one piece of text - paste() and paste0():
rownames(mv_cors) <- paste(tax_table(mv)[,6][1:50],(rownames(tax_table(mv)[1:50])))
colnames(mv_cors) <- paste(tax_table(mv)[,6][1:50],(rownames(tax_table(mv)[1:50])))


## Delete everything except the following, using setdiff() (difference between sets)
# via SO: https://stackoverflow.com/questions/2822532/how-can-i-neatly-clean-my-r-workspace-while-preserving-certain-objects
rm(list=setdiff(ls(), keepers))
# keepers can be your kept object or a list of objects to keep.. could also specify a subset of a list, like:
rm(list=setdiff(ls(), ls()[10:15]))


### TWO PROBLEMS: cant sort by two columns, and order keeps disappearing:
# kept getting bad orders of numeric data
# sort by column A, then by coulmn B:
dw_lefse_te <- dw_lefse_te[with(dw_lefse_te, order(Condition, -KW.RST)), ]
# however, still need to stop ggplot reordering:
# answer is to 'fix' the order by factorising:
dw_lefse_te$Taxa <- factor(dw_lefse_te$Taxa, levels(dw_lefse_te$Taxa))
#
## F I X   F A C T O R ?
# ## possibly set temps as a factor so keeps 'week' order? :
# mi.map$X.SampleID <- factor(mi.map$X.SampleID, levels = mi.map$X.SampleID)
## gavin simpson does it like this:
# theTable <- within(theTable, 
#                    Position <- factor(Position, 
#                                       levels=names(sort(table(Position), 
#                                                         decreasing=TRUE))))



## opposite to the %in% operator
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))


## grab significant species to a phylo object
## esp. U N I Q U E ( )
z.rm <- c(rownames(dw_d2_react_sig),rownames(dw_d2_tr1_sig),rownames(dw_d2_te_sig))
z.rm <- unique(z.rm)
z.rm2 <- which(rownames(otu_table(dw_slv)) %in% z.rm)
z.tax<-tax_table(dw_slv)[z.rm2,]
z.otu<-otu_table(dw_slv)[z.rm2,]
dw_slv_des <- merge_phyloseq(z.tax,z.otu,dw_slv_map)


## grab sums for taxa ranks w. phylogeny
z.rm <- tax_glom(dw_slv,'Phylum')
z.rm2 <- cbind(as(taxa_sums(z.rm),'numeric'),as(tax_table(z.rm))[,1],'data.frame')

## Clean up taxonomy table: set all  '' to NAs
## caffeine perdita

## !
## sort e.g. sample data by e.g. weeks  
sample_data(dw_slv)[with(sample_data(dw_slv),order(Week)),]


## Fix ACCESSION to REF-TAXONOMY
# R: need to negate the default quote reading, for unclear reasons (maybe blank ' "" ' header?)
z.con7 <- read.table('/media/cart/cart_sandbox/z.databases/silva128/taxonomy_97_16S/consensus_taxonomy_7_levels.txt',sep='\t',row.names=1, quote='')
library(tidyr)
z.con7 <- separate(data=z.con7, col=V2, into=c("Domain", "Phylum","Class","Order","Family","Genus","Species"), sep =';')
z.con7[,'Acc']<-rownames(z.con7)
z.recon <- as.data.frame(paste(z.con7$Domain, z.con7$Phylum, z.con7$Class, z.con7$Order, z.con7$Family, z.con7$Genus, z.con7$Species, z.con7$Acc, sep=';'))
row.names(z.recon) <- row.names(z.con7)
# make sure to put the object first, then the path!
write.table(z.recon,'consensus_taxonomy_7_levels_plus_Acc.txt',sep='\t')
# Term:
# nano: delete the column header; save.
# # remove the " " around the taxa
# sed 's/["]//g' '/media/cart/Dropbox/SilentGeno/R/R_DW/consensus_taxonomy_7_levels_plus_Acc.txt' > '/media/cart/Dropbox/SilentGeno/R/R_DW/consensus_taxonomy_7_levels_plus_Acc.txt'




# add minimal values for each column to each element in that column; effectively 'solving' negative values in a R-log/VST DF
# for z.rl<- assay(DESeq2::rlog(dataframe))
# get column minia, i.e. lowest transformed values for each sample  ::  simple additive pseudocount (not advisable!)
z.mini <- abs (apply (z.rl, 2, function(x) min(x) ))
z.mini
z.rl_fix <- as.data.frame( t (apply (z.rl, 1, function(x) x+z.mini )) )

# need to re-transpose as apply does some wierd transpose-and-cut
z.maxd<-t(apply(z.rl,1,function(x)x+z.mini))



# split columns in a DF (e.g. separate out taxa) - T I D Y R
library('tidyr')
TAX_DF<-(separate(data = BIG_DF, col = TAXCOL, into = c("Domain", "Phylum","Class","Order","Family","Genus","Species"), sep = "DELIM"))[,INDEX:INDEX]


# reverse or invert order of list / vector
rev(OBJECT)


##
## E D U A R D O
otu_plus_genes_higher = otu_plus_genes[(apply(otu_plus_genes[c(3:5)], 1, function(row) any(row > 1))), ]
## Generalise:
set_above1 = set [ (apply (set [c(3:5)], 1, function(row) any(row > 1)) ) , ]
z.rm5[apply(z.rm5,1,function(row) !any),]
## FOLLOWING THIS ADVICE< 
## to catch matching rows, perhaps instead use 'which' than 'any''any'
dw_ps[(which(dw_ps$OTU == 'sumden2')),]
##
## evaluate logic of 'x' in a range:
## check diff between '&' and '&&' in evaluation: 
## & : check each value in vector T or F
## && :check ALL values in vector are T / check ALL values in vector are F 
# as it happens, same in this case, as each matrix index is a single value:
sq.mat <- square matrix of values from -1 to 1 (e.g. correlation matrix)
range.TF <- apply(sq.mat,c(1,2),function(x) 0.5> x && x > -0.5)




##  G  G  P  L  O  T

#arbitrary legend placement
theme(..., legend.position =  c(.90, .95)) + 
  # change the theme of all plotting             <<  !  >>
  theme_classic() +
  # or by messing with the following
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey80") ) +
  # plotting by category requires being inside an aes() arg
  ggplot(DATAFRAME, aes(x=CATEGORY1, y=VALUES, color= CATEGORY2))
# horizontal line across axis
+ geom_hline(yintercept=0, color='COLOR#')
# change axis text size and orientation - lots of args!
+ theme(axis.text.x = element_text (angle=-90, hjust=0, vjust=0.5, size=14 ))
# title
+ labs(title='OVERALL RECRIMINIATION')
# another title for you
+ ggtitle("A TITLE FOR ALL SEASONS")
#rotate 90* or so..
+ coord_flip()
# change size by abundance (baseMean*2^(lfc)) in original DF
+ scale_size_continuous(name='SENSIBLE TITLE') + geom_point(alpha=0.6, aes(size=(baseMean*2^log2FoldChange))) 
# add centre points to balance alpha'd jitter
+ geom_point(color='COLOR',size=0.05)
# set colour legend to one column - works if remember to call on right aspect, e.g. fill/col/shape etc
+ guides(fill=guide_legend(ncol=3))
# alternatively, no legend box...
+ theme(legend.position="none")
# resize bars for equal space across TREATs:
+ facet_grid(~Reactor,scales='free_x',space='free_x')
# ... or across X and Y:
+ facet_grid(~Reactor,scales='free',space='free')
# ... subsetted facet of (A(B(C))) and Y:
+ facet_grid(A+B+C~Y)
# ...remove outlines (fill and outline set to same colour)...
+ geom_bar(aes(color=CATEG, fill=CATEG), stat ="identity", position="dodge")
# jitter the points to avoid overplotting, with transparency
+ geom_jitter(alpha=0.4, aes(size=sqrt(baseMean*2^log2FoldChange)))
# custom colour palette (see custom paletting below)
+ scale_color_manual(values = (z.col))
# custom FILLING colour, same palette
+ scale_fill_manual(values = (z.col))
# change BACKGROUND COLOUR, GID LINES
+ theme(panel.background = element_rect(fill='green', colour='red'))
theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
  theme(panel.grid.minor.y = element_line(colour='grey85', size = 0.1)) +                               # horiz lines 2
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  # change the facet grid (strip) boxes appearance
  theme(strip.background = element_rect(colour = NULL, fill = "white")) +
  # limit plotted colour values by setting what is displayed regardless of values plotted:
  + scale_colour_manual(values = cols, limits = c("4", "6", "8", "10"))
# manual assignment of NA values (while colours handled by colour vector:
+ scale_fill_manual(values=col.vector, na.value='grey50')                                                                  # custom-coloured bars 
# customise the shapes used in geom_point
+ scale_shape_manual(values = c(0, 5, 6, 15)) # NOTE: still need to assign a variable to 'shape' by

# the rotation of grid panel text                                               
+ theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270))
# increase width of geom_bar to avoid inconsistent gap betweenbars
+ geom_bar(width=1 ) # default apparently is 0.9  



## C O L O U R   P A L E T T E

# ## RColorBrewer
# brewer.pal(n, name)
# display.brewer.pal(n, name)
# display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
#                    colorblindFriendly=FALSE)
# brewer.pal.info


library(scales)
z.col<-hue_pal()(14)
z.col<-hue_pal()(length(get_taxa_unique(dw_full,'Genus')))
show_col(z.col)
# reverse the order of colours
show_col(rev(z.col))
z.col2<-sample(z.col)
show_col(z.col2)
plot_ordination(dw_full,ordinate(dw_full,'CCA',formula=~OLR),color='Phylum',shape='Treat',title='CCA of DW3up_Hell ~NH3',label='Timepoint',type='biplot') + scale_color_manual(values = z.col2)


## lots of colours
z.col <- colors(distinct=TRUE)[1:135] ; z.col <- c(z.col,colors(distinct=TRUE)[233:501])
z.col <- sample(z.col,85)   # pick 85 random colours
show_col(z.col)

# DW 12 colours (prob subset from the 42 colours below)
z.col12 <- c('#A51876','#E43FAD','#1E78D2','#094582','#117878','#3FE4E4',
             '#18A55E','#3FE491','#b8f9d3','#D2D21E','#D2781E','#E43F5B',
             'grey50')
show_col((z.col12))
z.col12 <- c(z.col12 )

# 39 colours:
z.col39 <- c("aquamarine4","goldenrod1","coral2","red","violetred","lightsteelblue3",
             "chocolate1","wheat1","purple","salmon4","mediumspringgreen","lightsalmon",
             "lightskyblue1","lightyellow4","darkolivegreen","tan1","lightpink1","tomato4","chocolate2",
             "oldlace","honeydew3","darkmagenta","pink1","brown1","darkgoldenrod4",
             "steelblue1","lightcoral","cadetblue4","greenyellow","azure","lightgoldenrod2",
             "darkkhaki","sienna1","orange2","olivedrab3","lightyellow","olivedrab2","rosybrown4",
             "peru")             
show_col(z.col39)

# 42 colours in cohesive pattern
z.col42 = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5",
            "#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA",
            "#98F0F0", "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811",
            "#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5","#784511","#A55E18",
            "#D2781E","#E4913F","#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C","#E43F5B",
            "#EA6C81","#F098A7")
show_col(z.col42)

# to 57 colours!
z.col57 <-c(z.col42,colours()[584:598])

z.col <- c("#781156","#114578","#117878","#18A55E","#A5A518","#784511","#781122","#F098A7","#185EA5","#3FE4E4")
