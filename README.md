## `<!> this document is not finished </!>`

This guide to metagenomic analysis is still being written (0ct '22). 
All (+/-)feedback is welcome: simply throw objects/comments directly at me, or [drop us a line at the related repo](https://github.com/handibles/climber/issues) .


# `CLI`, `M`icro`b`ial `E`cology, and `R`

Some basic steps in microbial ecology, focusing on the processing of `2ndGen` Illumina `fastq` data, into either `amplicon` (e.g. 16S) or `metagenomic` (e.g. shotgun) datasets, followed by ecology-based analysis of the communities and patterns we find in that data.


### Metagenomic data (i.e. shotgun)

  * <a href="https://handibles.github.io/climber/documents/shotgun_assembly.html">metagenomic shotgun assembly - cheese and/or chicken data </a> 

Metagenomic shotgun assembly covers the following steps:

  0. Setting up your analysis - `bash` and friends
  1. Checking your sequence data - `FastQC` & `MultiQC`
  2. Sequencing QC - filtering and trimming your sequences - `Trimmomatic`
  3. Sequencing QC - purifying your sequences - `BowTie2`
  4. Metagenomic Community profiling - `Kraken2` & `Bracken`


We also move through <a href="documents/data_to_R.html">importing output from `Kaiju` or `Kraken2+Bracken` into `R`<\a>. This awkward step is finished, but not pretty.

  * importing data into `R` - generating a count matrix, taxonomic table, and phyloseq object from metagenomic data


This metagenomic workflow is also present in simple, `raw code` (the two workflows will be unified once writing is finished, so not currently up to date with the above).

  * <a href="documents/shotgun_assembly_raw.html">`raw code only of the metagenomic shotgun assembly`</a> - as above, less explanation, **incomplete**


### Amplicon data (e.g. 16S)

Not yet started. The initial steps (setup, get data, QC) are very similar in most cases (remember to cut off your primers!), but are followed by a denoising step (`DADA2`) and optionally an attempt to predict the metabolic capabilities of the communities at hand (`PICRUSt2`).


### Microbial Ecology (and `R`)

Not yet started. This is the real magic, and we get to make _pictures_.


---
  <a href="documents/climber_todo.html">`</a>

<img src="vis/ucc.png" width="150" align="center" /> | <img src="vis/teag.png" width="150" align="center" /> | <img src="vis/v1group.png" width="150" align="center"/>
