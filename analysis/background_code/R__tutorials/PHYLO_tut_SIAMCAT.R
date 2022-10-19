
## ================================================= QUICK START ================================================= ##

library(SIAMCAT)

data("feat_crc_zeller", package="SIAMCAT")
data("meta_crc_zeller", package="SIAMCAT")

feat.crc.zeller[1:3, 1:3]
dim(feat.crc.zeller)
head(meta.crc.zeller)

# consider changing test on top example (gender?)
label.crc.zeller <- create.label(meta=meta.crc.zeller,
                                 label='Group', case='CRC')

siamcat <- siamcat(feat=feat.crc.zeller,
                   label=label.crc.zeller,
                   meta=meta.crc.zeller)
show(siamcat)

siamcat <- filter.features(siamcat,
                           filter.method = 'abundance',
                           cutoff = 0.001)


## ASSOCIATION TESTING

siamcat <- check.associations(
  siamcat,
  sort.by = 'fc',
  alpha = 0.05,
  mult.corr = "fdr",
  detect.lim = 10 ^-6,
  plot.type = "quantile.box",
  panels = c("fc", "prevalence", "auroc"))

## MODEL BUILDING

# NORMALISE
siamcat <- normalize.features(
  siamcat,
  norm.method = "log.unit",
  norm.param = list(
    log.n0 = 1e-06,
    n.p = 2,
    norm.margin = 1
  )
)

# CROSS-VALIDATE
siamcat <-  create.data.split(
  siamcat,
  num.folds = 5,
  num.resample = 2
)

# TRAIN
siamcat <- train.model(
  siamcat,
  method = "lasso"
)

model_type(siamcat)

models <- models(siamcat)
models[[1]]

# PREDICTIONS
siamcat <- make.predictions(siamcat)
pred_matrix <- pred_matrix(siamcat)
head(pred_matrix)



### MODEL EVAL AND INTERP.
siamcat <-  evaluate.predictions(siamcat)

# EVALUATION PLOT
model.evaluation.plot(siamcat)

# INTERPRETATION PLOT
model.interpretation.plot(
  siamcat,
  fn.plot = 'interpretation.pdf',
  consens.thres = 0.5,
  #norm.models = TRUE,
  limits = c(-3, 3),
  heatmap.type = 'zscore',
  prompt = FALSE
)

## ================================================= MODEL APPLICATION  ================================================= ##

## create (or recreate) the trained models on Zeller's French data 
siamcat.fr <- siamcat

## create siamcat obj from Yu's Chinese data
  fn.feat.cn  <-  'https://embl.de/download/zeller/CN-CRC/CN-CRC-N128_tax-ab-specI.tsv'
  fn.meta.cn  <-  'https://embl.de/download/zeller/CN-CRC/CN-CRC-N128_metadata.tsv'
  feat.cn  <- read.table(fn.feat.cn, sep='\t', quote="", check.names = FALSE)
  feat.cn.rel <- prop.table(as.matrix(feat.cn), 2)
  # metadata
  meta.cn  <- read.table(fn.meta.cn, sep='\t', quote="", check.names=FALSE, stringsAsFactors = FALSE)
  # SIAMCAT object
  siamcat.cn <- siamcat(feat=feat.cn.rel, meta=meta.cn, label='Group', case='CRC')
  siamcat.cn <- normalize.features(siamcat.cn,
                                   norm.param=norm_params(siamcat.fr),
                                   feature.type='original',
                                   verbose = 2)

## test models on Yu's Chinese data
  
## NORMALISE (taking exact method from french set)
  siamcat.cn <- normalize.features(siamcat.cn,
                                   norm.param=norm_params(siamcat.fr),
                                   feature.type='original',
                                   verbose = 2)
## PREDICT
siamcat.cn <- make.predictions(
  siamcat = siamcat.fr,
  siamcat.holdout = siamcat.cn,
  normalize.holdout = FALSE)

## EVALUATE
model.evaluation.plot('FR-CRC'=siamcat.fr,
                      'CN-CRC'=siamcat.cn,
                      colours=c('dimgrey', 'orange'))






