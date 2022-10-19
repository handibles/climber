
# for a feat_table j, metadata, and taxo, melt into columns

      ##  > >
      ##  > >  S O L U T I O N !!
      ##  > >
      
      jjj <- melt(cbind(j_crush, j[, jdh$order]), variable_name = "Sample")    # note fixing sort-by-order
      head(jjj)
      # take a variable: repeat each value ntaxa times
      a <- unlist(lapply(metadata[ jdh$order , "sample"], function(m) rep(m, nrow(taxo))))
      jjj <- cbind(jjj, a)
      head(jjj)
    

      
        ##  > >
        ##  > >  O T H E R   S O L U T I O N !!     :((((((((((((
        ##  > >
        
  
      ## SHRINK REDUNDANT CRUSHED, TAXA
      ## ERROR :: think this will just remove the first instance?...
        ## check length of uniques to be sure
          # # destroy tax, but remove dupe Gs
          # j_crush_reduced <- j_crush[!duplicated(j_crush[, "G"]) , ]
  
        #                        crush_tax                       feats
        jt <- cbind.data.frame(  j_crush_reduced, j_reduced[, order(m$sample)],
                               stringsAsFactors=FALSE )
        # View(jt)
        melt_jt <- melt(jt)
  
        # now join with sample data melted on variable (sample ID)
        melt_jtm <- (merge(melt_jt, m, by.x="variable", by.y="sample"))
        # View(melt_jtm)

      