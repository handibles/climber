# `CLI`, `M`icro`b`ial `E`cology, and `R`


Some basic steps in microbial ecology, focusing on the processing of `2ndGen` Illumina `fastq` data, into either `amplicon` (e.g. 16S) or `metagenomic` (e.g. shotgun) datasets, followed by ecology-based analysis of the communities and patterns we find in that data.

<img src="/vis/ucc.png" width="100" />

<img src="/vis/teag.png" width="100" />

<img src="/vis/v1group.png" width="100" />


## guides so far:

### Metagenomic data (i.e. shotgun)

Combining all the steps:

  * <a href="documents/shotgun_assembly.html">metagenomic shotgun assembly - cheese and/or chicken data</a>


The same process as above but separated into the different steps:

  0. <a href="/climber/documents/0.setup.html">Setting up your analysis</a>
  1. <a href="/climber/documents/1.checkdata.html">Checking your sequence data</a>
  2. <a href="/climber/documents/2.seqcleaning.html">Sequencing QC - filtering and trimming your sequences</a>
  3. <a href="/climber/documents/3.seqpurity.html">Sequencing QC - purifying your sequences</a>


### Amplicon data (e.g. 16S)

Not yet written. The initial steps (setup, get data, QC) are very similar in most cases (remember to cut off your primers!), but are followed by a denoising step (`DADA2`) and optionally an attempt to predict the metabolic capabilities of the communities at hand (`PICRUSt2`).


### Microbial Ecology (and `R`)

Not yet written. This is the real magic, and we get to make _pictures_.




  <a href="documents/climber_todo.html">`</a>
