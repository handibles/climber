
## igor bromhead


## bash!
  
  # sudo apt-get install openssl-dev        # necess for git2r, a knitr requirement
  # sudo apt-get install libmagick++-dev    # necess for magick, a phytools requirement
  # sudo apt-get install libglu1-mesa-dev   # for rgl, for clusterSim
  # sudo add-apt-repository -y ppa:cran/imagemagick
  # sudo apt install cmake                  # dependencies, you prob have this already
  # sudo apt-get install libmagick++-dev
  ## some is also for risk of rain...but thats another story

  # note : nonpareil uses R, but does not have a library.


## get the stuff we use for work
  
  install.packages("BiocManager")
  install.packages("devtools")
  
  
  # lapply it (deps for maaslin)
  library(BiocManager)   # biocLite deprecated
  deps <- c("ade4","agricolae","ALDEx2","ape","circlize","cluster",
            "clusterSim","ComplexHeatmap","compositions","cowplot",
            "dada2","DECIPHER","DESeq2","devtools","dynamicTreeCut",
            "FSA","gam","gamlss","gbm","ggplot2","ggpubr","ggvegan",
            "glmnet", "igraph","inlinedocs","knitr","logging",
            "MASS","networkD3","nlme","optparse","outliers","parallel",
            "penalized","phangorn","philr","phyloseq","propr","pscl",
            "RColorBrewer","reshape","robCompositions","robustbase",
            "scales","seqTools","shiny","shinythemes","tibble","vegan",
            "viridis","WGCNA","zCompositions", "nonpareil", "lme4",
            "abdiv", "rtk", "emmeans", "e1071")
  deps <- sort(unique(deps))
  lapply(deps, function(x) install(x, update = FALSE))


## get stuff from github:

  ## hmm... needs libcurl4
  
  # devtools::install_github("DeveloperName/PackageName")
  # "zachcp/phyloseq-tools",
  devs <- c("benjjneb/decontam",
            "gavinsimpson/ggvegan",
            # GMPR decommissioned?
             # "kassambara/ggpubr", # bioc
             "joey711/phyloseq",
             "liamrevell/phytools",
             "zachcp/phyloseq-tools")
  library('devtools')
  lapply(devs, function(x) devtools::install_github(x) )
  

##    < ! >    b u g s

  ## what was the fix for this again?...
  # installation path not writeable, unable to update packages: foreign, lattice, MASS, Matrix, mgcv, survival

  q('no')
