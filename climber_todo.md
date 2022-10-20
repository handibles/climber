# CLIMBER - CLI, Microbial Ecology, and R


Some basic steps in microbial ecology, focusing on the processing of `NextGen` illumina data, into either `amplicon` (e.g. 16S) or `metagenomic` (e.g. shotgun) datasets. 

If this appeals to you (and why wouldn't it), please dampen your squibs slightly - it's not written yet, and certianly barely ready to be used. Realistically, we'll need to separate things out into `recipes`, and have them cross-linked. This perhaps could be followed up with a `justthecode` flow for those who are keen to get thigns very, very wrong.

Enough of this; lets make like tracked changes and train.

```
      preamble
        ¬| see list of topics in analysis/*preamble*

      shotgun assembly
        ¬| assembly     : all the standard things to do

      amplicon assembly
        ¬| assembly     : all the standard things to do

```


## guides so far:

### Metagenomic Assembly

  * <a href="analysis/shotgun_assembly.html">metagenomic shotgun assembly - cheese and/or chicken data :chicken: :cheese: :DNA: </a>


## guides yet to start:

Should eventually cover:

### Metagenomic Assembly

  * <a href="">command-line intro and overview: `variables, for-loops, basic commands, general ideas etc.`</a>
  
  * <a href="">Quality checking: `FastQC and MultiQC`</a>
  * <a href="">Sequence Quality Control: `Trimmomatic` & friends</a>
  * <a href="">Sequence Purity: removing hosts and contaminants with `kneaddata` (& friends, if any)</a>
  
  * <a href="">Microbial Communities: `Kraken2+Bracken`, & friends</a>
  
  * <a href="">Metabolic Capacity: a general approach</a>
  * <a href="">Metabolic Capacity: The `HUMAnN` Pipeline</a>
  * <a href="">Metabolic Capacity: Anti-Microbial Resistance, Virulence, carbohydrate metabolism, and others</a>


### Microbial Community Analysis

  * if you thought the above list was bad...
