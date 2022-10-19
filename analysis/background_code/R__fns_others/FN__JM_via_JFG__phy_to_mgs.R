phy_to_mgs <- function (physeq, ...){
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  countData = round(as(otu_table(physeq), "matrix"), digits = 0)
  if (!is.null(sample_data(physeq, FALSE))) {
    ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))
  }
  else {
    ADF = NULL
  }
  if (!is.null(tax_table(physeq, FALSE))) {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
                                        data.frame(tax_table(physeq), stringsAsFactors = FALSE), row.names = taxa_names(physeq)))
  }
  else {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
                                        row.names = taxa_names(physeq)))
  }
  if (requireNamespace("metagenomeSeq")) {
    mrobj = metagenomeSeq::newMRexperiment(counts = countData, 
                                           phenoData = ADF, featureData = TDF, ...)
    if (sum(colSums(countData > 0) > 1) < ncol(countData)) {
      p = suppressMessages(metagenomeSeq::cumNormStat(mrobj))
    }
    else {
      p = suppressMessages(metagenomeSeq::cumNormStatFast(mrobj))
    }
    mrobj = metagenomeSeq::cumNorm(mrobj, p = p)
    return(mrobj)
  }
}
