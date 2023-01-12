
      rm(list=ls())
      set.seed(2205)


## 1 ======================================

    ## Kraken2+Bracken
      
      ## < ! >   make sur ethe path to the file matches your computer
  
      # output from kraken2-bracken, combined to one file
      krak2_brack_file="input/r0937/fhi__climber__nur_krakenStnd.15.3.10_abundances.tsv"    # abundance
      # output from single kraken2 report, converted to MPA
      krak2_mpa_file="input/r0937/fhi__climber__nur_krakenStnd.15.3.10_taxonomy.tsv"         # taxonomy
      
  
    ## kraken2 on its own
      
      just_krak <- "input/r0937/fhi__climber__nur_krakenMaxi15.3.10_kreportCombo_count.tsv"   # combined kraken
      just_krak_tax <- "input/r0937/fhi__climber__nur_krakenMaxi15.3.10_taxonomy.tsv"   # combined kraken taxonomies
      
          
    ## Kaiju
      
      # output from kaiju's kaiju2table output
      kaij_file="input/r0937/fhi__climber__nur_kaiju_summary_species.tsv"
    
    
## 2 ====================================
    
  ## a ====================================
      
      k2dat <- read.table(file = krak2_brack_file, header=TRUE, quote = "", sep = '\t')
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
    
      ## start with data from the abundance table names  --------------------------------
      
      tax_a <- data.frame(taxon=k2dat[ ,1], k2_id=k2dat[ ,2], row.names = rownames(k2dat) , stringsAsFactors = FALSE)
      tax_a[,1] <- gsub('\\/','', tax_a[,1])
      tax_a[,1] <- gsub('ALOs_','ALOs-', tax_a[,1])
      # tax_a[,1] <- gsub("'","", tax_a[,1])      # none?
      
      
      ## NOTE use quote="" to avoid nightmare of special characters like '"`/,| etc in names
      tax_b <- read.table( krak2_mpa_file, sep='\t', header=FALSE, fill = TRUE, stringsAsFactors = FALSE, quote="")
      head(tax_b) ; dim(tax_b)
      
      ## regularise taxonomic ranks (not all have K:P:C:O:F:G:S etc.)   -----------------
      
      ranks <- c("k__", "p__", "c__", "o__" , "f__", "g__", "s__")
      tax_b <- t(apply(tax_b, 1, function(aa){  # aa <- tax_b[1,]
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
      
      
      ## stitch  &  name   ----------------------------------
      
      tax_c <- merge(tax_b, tax_a, by.x="S", by.y="taxon")[, c(2:7,1,8)]
      head(tax_c) ; dim(tax_c)
      prenames_tax <- paste0("k2_", stringr::str_pad(tax_c[,"k2_id"], width=7, pad=0))
      rownames(tax_c) <- prenames_tax
      
      
      ## check all on the same page   -----------------------
      
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
      k2_tax <- tax_c[ pos_abund , ]
      
      mgdat <- data.frame("sample" = colnames(k2_counts),
                          "type" = gsub("(.).*", "\\1", colnames(k2_counts), perl = TRUE),
                          "seqdepth" = colSums(k2_counts) )
      # head(mgdat)
      # View(mgdat)
      
      ## check - do all match?   -----------------------
      
      dim(k2_counts)
      dim(k2_relab)
      dim(k2_tax)
      dim(mgdat)
      
      # look
      k2_counts[1:10, 1:10]
      
      # View(k2_counts)
      # hist( log10(k2_counts ))

            
  ## b ====================================
      
      # #   i n   b a s h   ================================================== ##
        # 
        # # presume you've made a load of kraken2 files:
        # ls -lsh $KROUT/*kraken_report
        # 
        # # be careful here that you don't accidentally include the $TEST file too! This would duplicate that sample
        # ~/bin/KrakenTools/combine_kreports.py -r $KROUT/*kraken_report -o $KROUT/krakenStnd_kreportCombo.tsv
        # 
        # # table of sorts
        # grep -E 'perc\stot_all' $KROUT/krakenStnd_kreportCombo.tsv | sed 's/#//' > $KROUT/krakenStnd_kreportCombo_count.tsv
        # grep -E '\sS1?\s' $KROUT/krakenStnd_kreportCombo.tsv >> $KROUT/krakenStnd_kreportCombo_count.tsv
      
      krak2 <- read.table(just_krak, header=TRUE, sep="\t", quote = "")
      rownames(krak2) <- paste0("k2_", stringr::str_pad(krak2$taxid, width=7, pad=0))
      head(krak2)
      # View(krak2)
      
      # choose to keep the abundances at "_lvl", rather than at "_all" 
      kfeat <- krak2[ , grep("\\d_lvl$", colnames(krak2), value= TRUE, perl=TRUE) ]
      head(kfeat)
      
      # fix some names, esp for the next step
      tax_x <- data.frame(taxon=krak2[ ,30], k2_id=krak2[ ,29], row.names = rownames(krak2) , stringsAsFactors = FALSE)
      tax_x[,1] <- gsub('\\/','', tax_x[,1])
      tax_x[,1] <- gsub('ALOs_','ALOs-', tax_x[,1])
      tax_x[,1] <- gsub('[\\[\\]]','', tax_x[,1], perl = TRUE)
      
      # make a clean genus and species column
      tax_x$G <- gsub("^\\s*(\\w*)\\s.*", "\\1", tax_x$taxon, perl=TRUE)
      tax_x$S <- gsub("^\\s*(\\w+)\\s(\\w+)\\s*.*", "\\2", tax_x$taxon, perl=TRUE)
      head(tax_x)
      # View(tax_x)
      
      # now gather abundances based on perfect matches to each "Genus species" combo
      kIDs_per_unique_binomials <- apply( unique(tax_x[,c("G", "S")]), 1, function(aa){ rownames(tax_x)[(tax_x$G == unlist(aa[1])) & (tax_x$S == unlist(aa[2]))] })
      k2only_counts <- t(sapply( kIDs_per_unique_binomials, function(aa){ colSums(kfeat[ aa , ]) }))
      
      ## note - here we impose a threshold of 10 reads total
      k2only_kIDs <- rownames(k2only_counts)[ rowSums(k2only_counts) > 10 ]
      
      ## make a taxonomy file: however, names dont match
      head(tax_a) ; dim(tax_a) #; View(tax_a)
      # tax_b is the mpa outout which reports all rows in the DB
      head(tax_b) ; dim(tax_b) #; View(tax_b)
      # however, maxikraken is a bigger DB!! 
      head(krak2) ; dim(krak2) #; View(krak2)
      # 140-144 matches
      sum( paste0("s__" , gsub(" ", "_", gsub("^\\s*", "", tax_x[ k2only_kIDs, "taxon"], perl=TRUE))) %in% k2_tax$S)
      sum( k2only_kIDs %in% rownames(k2_tax))
      
      # hands-on matching: 
      # in:
      paste0("s__" , gsub(" ", "_", gsub("^\\s*", "", tax_x[ k2only_kIDs, "taxon"], perl=TRUE)))[ paste0("s__" , gsub(" ", "_", gsub("^\\s*", "", tax_x[ k2only_kIDs, "taxon"], perl=TRUE))) %in% k2_tax$S ]
      # out:
      paste0("s__" , gsub(" ", "_", gsub("^\\s*", "", tax_x[ k2only_kIDs, "taxon"], perl=TRUE)))[ !(paste0("s__" , gsub(" ", "_", gsub("^\\s*", "", tax_x[ k2only_kIDs, "taxon"], perl=TRUE))) %in% k2_tax$S) ]
      
      ## match by genus? we have the species already
      unlist(
        lapply( tax_x[ k2only_kIDs, "taxon"], function(aa){   # aa <- tax_x[ k2only_kIDs, "taxon"][615] 
          bb <- strsplit(split = " ", gsub("^\\s*", "", aa, perl=TRUE))
          ifelse( grepl("virus", unlist(bb)), paste(bb, collapse = "_"), bb[[1]]) 
            })) 
            
      k2only_counts <- k2only_counts[ k2only_kIDs , ]
      k2only_tax <- tax_c[ k2only_kIDs , ]
      
      # check dimensions / orientation
      dim(k2only_counts)
      dim(k2only_tax)
      
      
  ## c ====================================
      
      # note the quote option, to stop strange names doing strange things
      kaidat <- read.table(file = kaij_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
      # remove 1+ header(s)
      kaidat <- kaidat[ !grepl("file", kaidat[,1]) , ]
      dim(kaidat)
      head(kaidat)
      # View(kaidat)
      
      ## fix sample names
      kaidat$file <- gsub(".*/(.*)_kaiju.*", "\\1", kaidat$file)
      kaidat$otu <- gsub( ".*;(.*);", "\\1", kaidat$taxon_name )
      kaidat$ID <- paste0( "kai__", stringr::str_pad( kaidat[, "taxon_id"], width = 7, pad = 0) )
      length(unique(kaidat$file))         # 20 files
      length(unique(kaidat$otu))          # 11.7K OTUs
      length(unique(kaidat$ID))           # 11.7K IDs
      kaidat$percent <- as.numeric(as.character(kaidat$percent))
      
      # should be very close to nsamples * 100
      sum(kaidat$percent, na.rm = TRUE)
      
      head(kaidat) # View(kaidat)
      str(kaidat)
      
      
  ## assess total output    ----------------------------
      
      kaidf <- as.data.frame(kaidat) 
      kaidf[ , "reads"] <- as.numeric(kaidf[ , "reads"])
      
      # again, counts and percent
      sum( (kaidf[ , "reads"]))
      sum( (kaidf[ , "percent"]), na.rm = TRUE)  # nsample*100%
      
      ## nreads unclassified / sent to bin
      kaidf$substrate <- gsub("([A-Z]).*", "\\1", kaidf$file)
      # ggplot(kaidf, aes(log10(reads), fill = substrate)) + geom_boxplot() + labs(title = "total N reads per substrate")
      
      
  ## congeal  -------------------------------------------
      
      # use the Kaiju ID, bot the genus name
      u_bio <- unique( kaidf$ID)
      u_sample <- unique( kaidf$file)
      
      ## count data
      # to each unique taxon found (u_bio), assign the number of counts found in each unique sample (u_sample)
      ka_counts <- do.call("rbind", lapply( u_bio, function(aa){   # aa <- u_bio[100]
        bb_df <- kaidf[ grep(aa, kaidf$ID) , ]
        cc_vec <- unlist(lapply( u_sample, function(aaa){ # aaa <- u_sample[9]
          as.numeric(bb_df[ match(aaa, bb_df$file) , "reads" ])
        }))
      } )
      )
      rownames(ka_counts) <- u_bio # [1:100]
      colnames(ka_counts) <- u_sample
      
      # only keep real numbers
      ka_counts[ is.na(ka_counts)] <- 0
      dim( ka_counts <- ka_counts[ , !apply( ka_counts, 2, function(aa){ all(aa == 0)}) ] )
      
      
      ## % data    
      # same, but make relative abudances using the values given by kaiju
      ka_relab <- do.call("rbind", lapply( u_bio, function(aa){   # aa <- u_bio[100]
        bb_df <- kaidf[ grep(aa, kaidf$ID) , ]
        cc_vec <- unlist(lapply( u_sample, function(aaa){ # aaa <- u_sample[9]
          as.numeric(bb_df[ match(aaa, bb_df$file) , "percent" ])
        }))
      } )
      )
      rownames(ka_relab) <- u_bio # [1:100]
      colnames(ka_relab) <- u_sample
      
      # only keep real numbers
      ka_relab[ is.na(ka_relab)] <- 0
      dim( ka_relab <- ka_relab[ , !apply( ka_relab, 2, function(aa){ all(aa == 0)}) ] )
      
      
      # look
      ka_counts[1:10,1:10]
      ka_relab[1:10,1:10]
      
      ka_tax <- unique( t(# do.call("rbind", 
        apply( kaidf[ , c("ID", "taxon_name")], 1, function(aa){   # aa <- kaidf[ 1, c("ID", "taxon_name" )]
          bb_vec <- aa[2]
          if( grepl("Viruses", bb_vec)){ 
            cc_vec <- c( rep( "Viruses", length(2:7)), aa[1])
          }else{
            cc_vec <- c( unlist(strsplit(unlist(bb_vec), ";"))[ 2:7], aa[1] )
          }
          cc_vec  }) ))
      
      
      # set names
      rownames(ka_tax) <- ka_tax[,7]
      colnames(ka_tax) <- c("p", "c", "o", "f", "g", "s", "ID")   # Jun'22 - remove "d" here, increase index of all names by 1
      head(ka_tax)
      dim(ka_tax <- (unique(ka_tax)))

      
## 3 ========================================================
      
      ## general metadata
      write.table(mgdat, "output/r0937__metadata__dat.tsv", sep = "\t")
      saveRDS( mgdat, "output/r0937__metadata__dat.RDS")
      
      
      ## for Kraken2+Bracken:
      
      ## tsv output
      write.table(k2_counts, "output/r0937__krakenStnd15.3.10__num.tsv", sep = "\t")
      write.table(k2_relab, "output/r0937__krakenStnd15.3.10__frac.tsv", sep = "\t")
      write.table(k2_tax, "output/r0937__krakenStnd15.3.10__tax.tsv", sep = "\t")
      
      ## RDS for R output - note samples are ROWS
      saveRDS( t(k2_counts), "output/r0937__krakenStnd15.3.10__num.RDS")
      saveRDS( t(k2_relab), "output/r0937__krakenStnd15.3.10__frac.RDS")
      saveRDS( k2_tax, "output/r0937__krakenStnd15.3.10__tax.RDS")
      
      
      ## for Kraken2-only:
      # huff
      dim(k2only_relab <- apply( k2only_counts, 2, function(aa) aa/sum(aa)))
      
      ## tsv output
      write.table(k2only_counts, "output/r0937__krakenMaxi15.3.10__num.tsv", sep = "\t")
      write.table(k2only_relab, "output/r0937__krakenMaxi15.3.10__frac.tsv", sep = "\t")
      write.table(k2only_tax, "output/r0937__krakenMaxi15.3.10__tax.tsv", sep = "\t")
      
      ## RDS for R output - note samples are ROWS
      saveRDS( t(k2only_counts), "output/r0937__krakenMaxi15.3.10__num.RDS")
      saveRDS( t(k2only_relab), "output/r0937__krakenMaxi15.3.10__frac.RDS")
      saveRDS( k2only_tax, "output/r0937__krakenMaxi15.3.10__tax.RDS")
      
      
      ## for Kaiju:
      
      ## tsv output
      write.table(ka_counts, "output/r0937__kaijuStnd__num.tsv", sep = "\t")
      write.table(ka_relab, "output/r0937__kaijuStnd__frac.tsv", sep = "\t")
      write.table(ka_tax, "output/r0937__kaijuStnd__tax.tsv", sep = "\t")
      
      ## RDS for R output - note samples are ROWS
      saveRDS( t(ka_counts), "output/r0937__kaijuStnd__num.RDS")
      saveRDS( t(ka_relab), "output/r0937__kaijuStnd__frac.RDS")
      saveRDS( ka_tax, "output/r0937__kaijuStnd__tax.RDS")
      
      
## 4  =========================================================
      
      ## note: we use the count data (making % is easy in R, but 
      ##  - we would have no way of knowing the original counts)
      
      ## note also: we did not extract any phylogenetic tree from 
      ##  - the shotgun data, so no tree is included below
      
      library(phyloseq)
      
  ## Kraken2/Bracken   ---
      
      # look at the first ten rows, ten columns:
      k2_counts[1:10, 1:10]
      k2_relab[1:10, 1:10]
      k2_tax[1:10, 1:5]
      mgdat[1:10,  ]
      
      # note: we use the count data
      k2_otu <- otu_table(k2_counts, taxa_are_rows = TRUE)
      # convert from df to matrix to tax_table to avoid losing rownames...
      k2_tax2 <- tax_table(as.matrix(k2_tax))
      
      (k2_phyloseq_dat <- phyloseq(
        k2_otu,
        k2_tax2,
        # note no phylogenetic tree
        sample_data(mgdat)
      ))
      
      saveRDS( k2_phyloseq_dat, "output/r0937__KrakenStnd15.3.10__phylo.RDS")
      # check: save and read back in. RDS is a handy way to save/load data
      k2_phylo <- readRDS("output/r0937__KrakenStnd15.3.10__phylo.RDS")
      k2_phylo
      
      
  ## Kraken2 ONLY   ---
      
      # look at the first ten rows, ten columns:
      k2only_counts[1:10, 1:10]
      k2only_relab[1:10, 1:10]
      k2only_tax[1:10, 1:5]
      mgdat[1:10,  ]
      
      # note: we use the count data
      k2only_otu <- otu_table(k2only_counts, taxa_are_rows = TRUE)
      # convert from df to matrix to tax_table to avoid losing rownames...
      k2only_tax2 <- tax_table(as.matrix(k2only_tax))
      
      (k2only_phyloseq_dat <- phyloseq(
        k2only_otu,
        k2only_tax2,
        # note no phylogenetic tree
        sample_data(mgdat)
      ))
      
      saveRDS( k2only_phyloseq_dat, "output/r0937__KrakenMaxi15.3.10__phylo.RDS")
      # check: save and read back in. RDS is a handy way to save/load data
      k2only_phylo <- readRDS("output/r0937__KrakenMaxi15.3.10__phylo.RDS")
      k2only_phylo
      
      
  ## Kaiju   ----------
      
      # look at the first ten rows, ten columns:
      ka_counts[1:10, 1:10]
      ka_relab[1:10, 1:10]
      ka_tax[1:10, 1:5]
      mgdat[1:10, ]
      
      # note: we use the count data..
      # sample names need to match - here we solve an issue between "-" and "."
      colnames(ka_counts) <- gsub("-", ".", colnames(ka_counts))
      dim(kai_otu <- otu_table(ka_counts, taxa_are_rows = TRUE))
      # minimum modification of sample names
      colnames(kai_otu) <- paste0("X", colnames(kai_otu))
      dim(kai_tax2 <- tax_table(ka_tax))
      
      (kai_phyloseq_dat <- phyloseq(kai_otu, kai_tax2, sample_data(mgdat)))
      
      saveRDS( kai_phyloseq_dat, "output/r0937__KaijuStnd__phylo.RDS")
      # demo - save out and read back in again
      kai_phylo <- readRDS("output/r0937__KaijuStnd__phylo.RDS")
      kai_phylo