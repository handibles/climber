# some propr seaweed 
library('phyloseq')
library('propr')
library('zCompositions')
PHYLO

# # overhead tut stuff
# set.seed(12345)
# N <- 100
# X <- data.frame(a=(1:N), b=(1:N) * rnorm(N, 10, 0.1),
#                 c=(N:1), d=(N:1) * rnorm(N, 10, 1.0))


## transform first, then subset (although probably robust to such things)

# from base of PHYLO_commmcomp
d.1 <- t(data.frame(otu_table(PHYLO)))
dim(d.1)
# choice of 0-wrangling here important :: GBM v. CZM - joey711 argues GBM not applicable to micro (i.e. no real zeroes)
x <- cmultRepl(d.1, method='CZM', output = 'p-counts')   #  SAMPLES AS ROWS


# x <- as.data.frame(t(otu_table(PHYLO)))
# #x[1:2,1:2]

phi <- propr(x, metric = "phi", symmetrize = TRUE)
rho <- propr(x, metric = "rho", ivar = 0)
phs <- propr(x, metric = "phs", ivar = 0)

# ?
updateCutoffs(rho, cutoff =  seq(.05, .95, .3))


# # dont, i think, do this, as not compo-subcoherent
# feb <- x[ , grep('FEB', colnames(x)) ]

## these methods are alternative methods; create an index or copy on modify
  # index-on-fly (indexed subset demarked within objects matrix slot) versus... 
  rho99 <- rho[">", .99999995]
  simplify(rho99)  # preferable to subset, see vignette
  
  # copy-on-modify
  # more useful for a clustering subset
  rhoab <- subset(rho, select = c("a", "b"))
  rhoab@matrix

## ===

# vis
plot(rho99)


