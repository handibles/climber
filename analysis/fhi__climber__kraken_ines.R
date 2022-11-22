      
      
      mfile="input/fhi_redch_krakenStnd_abundances.tsv"
      
      
      jj <- read.table(file = mfile, header=TRUE, quote = "", sep = '\t')
      jj[1:10, 1:10]
      
      prenames_jj <- paste0("k2_", stringr::str_pad(jj[,2], width=7, pad=0))
      ## non-unique names -  assume only one dupe of each
      prenames_jj[ duplicated(jj[,2])] <- gsub('k2_', 'k2-A_', prenames_jj[ duplicated(jj[,2])])
      rownames(jj) <- prenames_jj
      
      # counts
      j <- jj[ , grep('_num', colnames(jj))]
      colnames(j) <- gsub( ".bracken_num", "", colnames(j) )
      
      # %
      k <- jj[ , grep('_frac', colnames(jj))]
      colnames(k) <- gsub( ".bracken_frac", "", colnames(k) )
      
      
      ##   T A X A  ::   K r a k e n   O u t p u t
      
      tax_a <- data.frame(taxon=jj[ ,1], k2_id=jj[ ,2], row.names = rownames(jj) , stringsAsFactors = FALSE)
      tax_a[,1] <- gsub('\\/','', tax_a[,1])
      tax_a[,1] <- gsub('ALOs_','ALOs-', tax_a[,1])
      # tax_a[,1] <- gsub("'","", tax_a[,1])      # none?
      
      ## NOTE escape quotes to stop the world ending
      tax_b <- read.table('input/fhi_redch_krakenStnd_taxonomy.tsv', sep='\t', header=FALSE, fill = TRUE, stringsAsFactors = FALSE, quote="")
      head(tax_b)
      
      ## A R R E T E !! 
      # issue is that ranks are not set, so a taxon missing e.g. order has all ranks shifted up by one, despite having the correct preceding term, "o__"
      ranks <- c("k__", "p__", "c__", "o__" , "f__", "g__", "s__")
      tax_b <- t(apply(tax_b, 1, function(aa){  # aa <- tax_b[1,]
        sapply(ranks, function(aaa){ # aaa <- "c__"
          bbb <- grep(aaa, unlist(aa), value=TRUE)
          ifelse( length(bbb) == 0, "unkn.", bbb)
        })
      }))
      dim(tax_b) ; str(tax_b)
      
      ## eukarya have loads of ranks, so many empty sapces to be filled, e.g.
      # tax_b[4000:4200,]
      
      # more issues due to names: doubled ranks:
      prior_index <- (unlist( lapply(1:nrow(tax_b), function(aa){ if( sum( grepl("unkn.", tax_b[aa,])) == 6){aa} }) ) - 1)
      # # compare:
      # tax_b[ prior_index , ]
      # tax_b[ prior_index+1 , ]
      # amend, by doubling up (will remove duplicates later)
      tax_b[ prior_index , "s__"] <- tax_b[ (prior_index+1) , "s__"]
      # row is only useful if has a species: cannot have a species, and hit 6, And be worth keeping
      tax_b <- tax_b[ !apply(tax_b, 1, function(aa) sum(grepl("unkn.", unlist(aa)) ) == 6) , ]
      
      
      colnames(tax_b) <- c("D", "P", "C", "O", "F", "G", "S")
      
      ## correct the double rank in tax_b  -  not a problem in K2_STND ?
      # x <- c("d", "p", "c", "o" , "f", "g", "s")
      # for(a in 1:7){    # a <-3
      #   tax_b[,a] <- gsub(paste0(x[a],"_",x[a],"__"), paste0(x[a],"__"), tax_b[,a])       # YES FINE I USED A FORLOOP, JESUS CALM DOWN
      # }
      
      tax_a[ , "taxon"] <- gsub(" ", "_", paste0("s__", tax_a[ , "taxon"]))
      head(tax_a) ; head(tax_a[ , "taxon"]) ;  dim(tax_a)
      head(tax_b) ;   head(tax_b[ , "S"]) ;  dim(tax_b)
      
      
      ## stitch  &  name
      tax_c <- merge(tax_b, tax_a, by.x="S", by.y="taxon")[, c(2:7,1,8)]
      head(tax_c) ; dim(tax_c)
      prenames_tax <- paste0("k2_", stringr::str_pad(tax_c[,"k2_id"], width=7, pad=0))
      rownames(tax_c) <- prenames_tax
      
      
      ##   U N I F Y 
      
      # check all on the same page
      all(rownames(jj) == rownames(j) & rownames(jj) == rownames(k))
      all(rownames(jj) %in% rownames(tax_c))
      
      kraken_ids <- sort(rownames(jj))
      jj <- jj [kraken_ids , ]
      j <- j [kraken_ids , ]
      k <- k [kraken_ids , ]
      tax_c <- tax_c[kraken_ids , ]
      
      
      ## remove if no values (kraken2 set to give all taxa)
      pos_abund <- rownames(k)[rowSums(k) > 0]
      k <- k[ pos_abund , ]
      j <- j[ pos_abund , ]
      tax <- tax_c[ pos_abund , ]
      
      
      ##   I N V E N T   M E T A D A T A
      mgdat <- data.frame("sample" = colnames(j),
                         "type" = gsub("(.).*", "\\1", colnames(j), perl = TRUE),
                         "seqdepth" = colSums(j) )
      
      
      ## check - all match?
      dim(j)
      dim(k)
      dim(tax)
      dim(mgdat)
      
      j[1:10, 1:10]
      
      
      ## tsv output
      write.file(j, "output/fhi__redch__krakenStnd-20__num..tsv")
      write.file(k, "output/fhi__redch__krakenStnd-20__frac..tsv")
      write.file(tax, "output/fhi__redch__krakenStnd-20__tax..tsv")
      write.file(mgdat, "output/fhi__redch__krakenStnd-20__dat..tsv")
      
      # RDS for R output - note samples are ROWS
      saveRDS( t(j), "output/fhi__redch__krakenStnd-20__num.RDS")
      saveRDS( t(k), "output/fhi__redch__krakenStnd-20__frac.RDS")
      saveRDS( tax, "output/fhi__redch__krakenStnd-20__tax.RDS")
      saveRDS( mgdat, "output/fhi__redch__krakenStnd-20__dat.RDS")

      
## phyloseq, see http://joey711.github.io/phyloseq/import-data.html for the most fun
      
      # note: we use the count data..
      (phyloseq_dat <- phyloseq(j, tax, mgdat))
      # note we lack a tree - we did not extract phylogenetic tree from the shotgun data
      saveRDS( phyloseq_dat, "output/fhi__redch__krakenStnd-20__phylo.RDS")
      