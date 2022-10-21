# `CLI`, `M`icro`b`ial `E`cology, and `R`


Some basic steps in microbial ecology, focusing on the processing of `2ndGen` Illumina `fastq` data, into either `amplicon` (e.g. 16S) or `metagenomic` (e.g. shotgun) datasets, followed by ecology-based analysis of the communities and patterns we find in that data.


## guides so far:

### Metagenomic data (i.e. shotgun)

Combining all the steps:

  * <a href="documents/shotgun_assembly.html">metagenomic shotgun assembly - cheese and/or chicken data</a>


The same process as above but separated into the different steps:

  1. <a href="/documents/0.setup.html">0 - Setting up your analysis</a>
  2. <a href="/documents/1.checkdata.html">1 - checking your sequence data</a>
  3. <a href="/documents/2.seqcleaning.html">2 - sequencing QC - filtering and trimming your sequences</a>
  4. <a href="/documents/3.seqpurity.html">3 - sequencing QC - purifying your sequences</a>


### Amplicon data (e.g. 16S)

Not yet written. The initial steps (setup, get data, QC) are very similar in most cases (remember to cut off your primers!), but are followed by a denoising step (`DADA2`) and optionally an attempt to predict the metabolic capabilities of the communities at hand (`PICRUSt2`).


### Microbial Ecology (and `R`)

Not yet written. This is the real magic, and we ge to make _pictures_.



  <a href="documents/climber_todo.html">`</a>
