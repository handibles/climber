
  
## >>> you can if you like, and it will emphasise rarer contributions, but there needs to be a clear rationale
## >>> and You don't understand enough to provide that.

## testing sqrt values in Bray curtis

  source("analysis/ms__som__chunk__spreadsheeting_env.R")
  library("vegan")
  library("ggplot2")

  # make sqrt
  dim(mgfeat_sqrt <- t(apply(mgfeat, 1, sqrt)))
  
  mgfeat[1:10, 1:10]
  mgfeat_ra[1:10, 1:10]
  mgfeat_sqrt[1:10, 1:10]

  hist( rowSums(mgfeat))
  hist( rowSums(mgfeat_ra))
  hist( rowSums(mgfeat_sqrt)) # nto a stndardisation, so no need to apply()
  
  bc <- vegdist(mgfeat, method = "bray")
  bc_ra <- vegdist(mgfeat_ra, method = "bray")
  bc_sq <- vegdist(mgfeat_sqrt, method = "bray")
  
  pcoa_bc <- wcmdscale(bc)
  pcoa_bc_ra <- wcmdscale(bc_ra)
  pcoa_bc_sq <- wcmdscale(bc_sq)
  
  
  ggpubr::ggarrange(
    ord_plot(vegout = pcoa_bc, 
             pdatf = mgdat, featdf = mgfeat) + coord_cartesian(),
    ord_plot(vegout = pcoa_bc_ra, 
             pdatf = mgdat, featdf = mgfeat_ra) + coord_cartesian(),
    ord_plot(vegout = pcoa_bc_sq, 
             pdatf = mgdat, featdf = mgfeat) + coord_cartesian(),
    nrow = 2, ncol = 2, common.legend =TRUE
    )
  
  ##  
  ##  
  ##  
  # https://ordnews.colostate.narkive.com/lMWF502c/1593-log-sqrt-and-other-transformation-with-bray-curtis-dissimilarity
  ##  
  ##  
  ##  

    par(mfrow = c(1,2))  
  # slope should be high (density = abundance per sample?)
    a_v <- cbind(
      (scale(apply(mgfeat, 2, mean))),
      (scale(apply(mgfeat, 2, var)))
    )
    rownames(a_v) <- mgtax[ rownames(a_v), 6]
    colnames(a_v) <- c("mean","var")
    #
    lm(var ~ mean, data.frame(a_v))
    cor(a_v, method = "sp")
    plot(a_v, type = "n", main = "scaled stats on raw ab",
    sub = "Logically, when relative densities among taxa are high the emerging\n
    relationship between averages and variances will be clearly linear...",
    xlab = NA)
    text(a_v, labels = rownames(a_v), col = "red", cex = 0.9)

    
  # applying a transform should reduce the slope - it does not.
    a_vs <- cbind(
      scale(apply(sqrt(mgfeat), 2, mean)),
      scale(apply(sqrt(mgfeat), 2, var))
    )
    rownames(a_vs) <- mgtax[ rownames(a_vs), 6]
    colnames(a_vs) <- c("mean","var")
    #
    lm(var ~ mean, data.frame(a_vs))
    cor(a_vs, method = "sp")
    plot(a_vs, type = "n", main = "scaled stats on sqrt",
    sub = "...After applying a transformation and plotting averages versus\n
    variances for the transformed data, the slope will be reduced",
         xlab = NA)
    text(a_vs, labels = rownames(a_vs), col = "red", cex = 0.9)
  
    ## sama sama    
    ggpubr::ggarrange(
        ggplot( as.data.frame(a_v), aes(x = mean, y = var)) +
          geom_point() +
          annotate("text", label = paste0("slope = ",(coef(lm(var ~ mean, data.frame(a_v))))[[2]] ), x = 5, y = 15) +
          geom_smooth(method ="lm"),
        ggplot( as.data.frame(a_vs), aes(x = mean, y = var)) +
          geom_point() +
          annotate("text", label = paste0("slope = ",(coef(lm(var ~ mean, data.frame(a_vs))))[[2]] ), x = 5, y = 15) +
          geom_smooth(method ="lm")
    )

    print(paste(    
    "Bray-Curtis values will be dominated by dominant species if no transformation is
    applied to the raw data: whether that is what we are after is a matter of our
    research question." ))
        
