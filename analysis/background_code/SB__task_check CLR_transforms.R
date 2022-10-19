
# checking CLR stuff  -  must be cmult replaced first

    ## better place for this     
    x <- matrix(c(1,1,1,
        1,1,1,
        1,1,1,
        0,1,1), nrow=4)
    
    ## samples are ROWS for ALL
    
    xx <- cmultRepl(x)      
    
    xxx_j <- t(apply(xx, 1, function(x){log(x) - mean(log(x))}))
    
    xxx_s <-  sweep(log(xx), 1, rowMeans(log(xx)), "-")
    
    xxx_c <- compositions::clr(xx)
    
    all(     # three checks is tautology but
      xxx == xxx_s,
      xxx_s == xxx_c,
      xxx == xxx_c )

