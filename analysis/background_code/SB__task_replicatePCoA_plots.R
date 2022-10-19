## be able to separate out the different groups etc
DATA <- 
  GROUP <- "control"  
CATEG <- "treat"


too_many_pcoas <- function(DATA){
  
  ggplot( DATA,aes(x = PCoA1, y = PCoA2, shape=schedule)) +   #, color=visit
    stat_ellipse(aes(x = PCoA1,y =PCoA2, fill=visit, colour=visit, lty=schedule), geom="polygon", size=1.2, level=0.8 , alpha=0.0, weight=0.8) +   #, colour=event, colour=visit
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "grey80") ) +
    scale_fill_manual("visit",values=col_var) +
    scale_color_manual("visit",values=col_var ) +
    scale_shape_manual("schedule",values=c(21, 22, 24)) +
    scale_alpha_manual("schedule",values=c(1, 0.3, 0.8)) +
    scale_linetype_manual("schedule",values=c(3,1,2)) +
    #   geom_line(aes(group = subject), alpha=0.15, size=0.3, colour='grey28') +
    
    < !! >
  ## difficult here to set values for subsets prgrammatically
    < !! >

    geom_point(aes(fill=factor(DATA$visit), shape=factor(DATA$schedule), alpha=factor(DATA$schedule)), color="grey28", size = 3, stroke=0.8) +    # color=visit,   
    xlab(PCo1) +
    ylab(PCo2) + border(color = 'grey35') +
    theme(legend.position='left', plot.margin = margin(2, 2, 0, 0, 'cm'))
  
}


x <- (bc_df[ bc_df$treat == "Test" , ] )
too_many_pcoas(x)

too_many_pcoas(bc_df)
too_many_pcoas(bc_df[bc_df$treat != "Test" , ])
too_many_pcoas(bc_df[bc_df$treat != "Control" , ])

