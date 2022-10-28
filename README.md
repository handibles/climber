## `<!> this document is not finished </!>`

This guide to metagenomic analysis is still being written (0ct '22). 
All (+/-)feedback is welcome: simply throw objects/comments directly at me, or [drop us a line at the related repo](https://github.com/handibles/climber/issues) .


# `CLI`, `M`icro`b`ial `E`cology, and `R`

Some basic steps in microbial ecology, focusing on the processing of `2ndGen` Illumina `fastq` data, into either `amplicon` (e.g. 16S) or `metagenomic` (e.g. shotgun) datasets, followed by ecology-based analysis of the communities and patterns we find in that data.


## guides so far:

### Metagenomic data (i.e. shotgun)

  * <a href="https://handibles.github.io/climber/documents/shotgun_assembly.html">metagenomic shotgun assembly - cheese and/or chicken data </a> 

Metagenomic shotgun assembly covers the following steps:

  0. Setting up your analysis
  1. Checking your sequence data
  2. Sequencing QC - filtering and trimming your sequences
  3. Sequencing QC - purifying your sequences - **not yet finished**
  4. Kraken2 (**not yet started**)
  5. Bracken (**not yet started**)

A similar workflow is also present in simple, `raw code` (the two workflows will be unified once writing is finished).

  * <a href="documents/shotgun_assembly_raw.html">`raw code only of the metagenomic shotgun assembly`</a> - as above, but with less explanation 



### Amplicon data (e.g. 16S)

Not yet started. The initial steps (setup, get data, QC) are very similar in most cases (remember to cut off your primers!), but are followed by a denoising step (`DADA2`) and optionally an attempt to predict the metabolic capabilities of the communities at hand (`PICRUSt2`).


### Microbial Ecology (and `R`)

Not yet started. This is the real magic, and we get to make _pictures_.


---
  <a href="documents/climber_todo.html">`</a>

<img src="vis/ucc.png" width="150" align="center" /> | <img src="vis/teag.png" width="150" align="center" /> | <img src="vis/v1group.png" width="150" align="center"/>
