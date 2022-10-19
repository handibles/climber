
# via phyloseq's UniFrac handling
## PCoA

uf.pcoa <- ordinate(eoe, 'PCoA', 'unifrac', weighted=TRUE)  # note double call for UF, and weighting arg

## ================================================================================================

## labels: PCoA has pre-calcd this for you (relative eigen)
PCo1 <- paste0('PCoA1: ',round(uf.pcoa$values[[2]][1]*100, 0),'%')
PCo2 <- paste0('PCoA2: ',round(uf.pcoa$values[[2]][2]*100, 0),'%')
PCo3 <- paste0('PCoA3: ',round(uf.pcoa$values[[2]][3]*100, 0),'%')
PCo4 <- paste0('PCoA4: ',round(uf.pcoa$values[[2]][4]*100, 0),'%')
PCo5 <- paste0('PCoA5: ',round(uf.pcoa$values[[2]][5]*100, 0),'%')
PCo6 <- paste0('PCoA6: ',round(uf.pcoa$values[[2]][6]*100, 0),'%')
PCo7 <- paste0('PCoA7: ',round(uf.pcoa$values[[2]][7]*100, 0),'%')
PCo8 <- paste0('PCoA8: ',round(uf.pcoa$values[[2]][8]*100, 0),'%')
PCo9 <- paste0('PCoA8: ',round(uf.pcoa$values[[2]][9]*100, 0),'%')
PCo10 <- paste0('PCoA10: ',round(uf.pcoa$values[[2]][10]*100, 0),'%')

eoe_rich <- estimate_richness(eoe)   #duplicate
uf.df = data.frame(data.frame(sample_data(eoe), stringsAsFactors = FALSE), eoe_rich,
                   "PCoA1" = uf.pcoa$vectors[,1]*-1,
                   "PCoA2" = uf.pcoa$vectors[,2]*-1,
                   "PCoA3" = uf.pcoa$vectors[,3]*-1,
                   "PCoA4" = uf.pcoa$vectors[,4]*-1,
                   "PCoA5" = uf.pcoa$vectors[,5]*-1,
                   "PCoA6" = uf.pcoa$vectors[,6]*-1,
                   "PCoA7" = uf.pcoa$vectors[,7]*-1,
                   "PCoA8" = uf.pcoa$vectors[,8]*-1,
                   "PCoA9" = uf.pcoa$vectors[,9]*-1,
                   "PCoA10" = uf.pcoa$vectors[,10]*-1,
                   "PCoA11" = uf.pcoa$vectors[,11])

## ================================================================================================

# a-div  -  violin plot, grouped by col of choice
ggplot(uf.df,aes(x=Description,y=uf.df$Shannon,fill=Condition)) + 
  geom_violin() +
  scale_fill_manual("Description",values=c(e.col)) +
  geom_dotplot(binaxis="y",stackdir = "center",dotsize=1,binwidth = 0.03, fill="black") + # take black out of the aes bracket!
  theme_classic()

uf_pcoa_boxp = cbind(data.frame(sample_data(eoe)$Condition), data.frame(sample_data(eoe)$Treatment), data.frame(sample_data(eoe)$Description), uf.pcoa$vectors)
uf_pcoa_boxp <- melt(uf_pcoa_boxp)
colnames(uf_pcoa_boxp) <- c('Condition', 'Treatment' , 'Description','PC','value')

# note low number of axes
ggboxplot(uf_pcoa_boxp[1:(0.4*nrow(uf_pcoa_boxp)),], 'PC', "value", fill = 'Description', size=0.3, palette = e.desc, width = 0.5 ) +
  theme(axis.text.x = element_text (angle=-90, hjust=0, vjust=0.5, size=11 )) +   #, legend.position =  c(.90, .95)
  labs(title='PCoA eigenvalues (20 of 64) of treatment*pre/post')

## ================================================================================================

gg.pcoa <- ggplot(uf.df,aes(x = PCoA1, y = PCoA2, color=Description, shape=Description)) +
  stat_ellipse(aes(x = PCoA1,y =PCoA2, fill=Description), geom="polygon" , level=0.8 , alpha=0.2) +   #, linetype=Description
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey80") ) +
  scale_fill_manual("Description", values=e.desc) +
  scale_color_manual("Description", values=e.desc ) +
  scale_shape_manual("Description", values=c(21,24, 22, 22, 22)) +
  geom_line(aes(group = ID), alpha=0.65, size=0.3, colour='grey28') +
  geom_point(color='grey28', aes(fill=factor(uf.df$Description)), size = 4) +
  xlab(PCo1) +
  ylab(PCo2) + border(color = 'grey35') +
  theme(legend.position='left', plot.margin = margin(2, 2, 0, 0, 'cm'))

# put boxplots inside plot, not alongside.  
xbp <- ggboxplot(uf.df, "Description", "PCoA1", fill = "Description", size=0.3, palette = e.desc, width = 0.5 ) +
  rotate() +  theme_transparent() +  theme(legend.position='none')

ybp <- ggboxplot(uf.df, "Description", "PCoA2", fill = "Description", size=0.3, palette = e.desc, width = 0.5 ) +
  theme_transparent() +  theme(legend.position='none')

xbp_grob <- ggplotGrob(xbp) ; ybp_grob <- ggplotGrob(ybp)

# add spacer, insert
gg.pcoa + annotation_custom(grob = xbp_grob, 
                            xmin = -0.5, xmax = 0.5, 
                            ymin = 0.35, ymax = 0.48) +
  annotation_custom(grob = ybp_grob,
                    xmin = 0.48, xmax = 0.62, 
                    ymin = -0.4, ymax = 0.4)

## ================================================================================================
