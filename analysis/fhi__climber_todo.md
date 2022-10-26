---
title: 'climber - tasks outstanding'
output:
  html_document:
    output_file: "climber_todo.html"
    toc: FALSE
knit: (function(inputFile, encoding) { rmarkdown::render(inputFile, encoding=encoding, output_file='../documents/climber_todo.html') })
---

## `<!> this document is not finished </!>`

This guide to metagenomic analysis is still being written (0ct '22). 
All (+/-)feedback is welcome: simply throw objects directly at me, or drop us a line [at the related repo](https://github.com/handibles/climber/issues) .


#  CLI, Microbial Ecology, and R

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
  

#### Picassan Refernces

  * [microbiome helper](https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(Humann2)) - Langille et al. 2017-2021

  
#### Presentation

Probably would do the best / look the best running on jekyll wiht the readthedocs theme (nice in teal, no?)
  
  
#### About

Written 2022 with the help of Jamie FitzGerald, Ines Calvete, Nuria Pena-Vidal, and all or most of those at UCC's [Claesson Lab](http://github.com/ClaessonLabUCC), and Teagasc's [Vision1](https://www.teagasc.ie/news--events/daily/other/teagascs-vision-1-lab-wins-at-the-irish-laboratory-awards.php) group.

---


#### a quick note about `glomming` and `sed`

Lets look at the `sed` ("`s`tream `ed`itor") command again! `r emo::ji("smile")` It can be used in lots of ways, but one of the most popular is to use [regular expressions (`regex`)](https://robinwinslow.uk/regex), a set of categories that allow matching and searching for different characters. We won't go into all the details here, but it's worth explaining why we do this/how this part works. 

The terminal does something similar when we use `*`, which _in the terminal_ means "as many of any character as possible". This is called [**`glomming`**](https://linuxhint.com/bash_globbing_tutorial/) The command `ls ~/data/ref/bt2_indices/*fna.gz` will normally list all files in `ls ~/data/ref/bt2_indices/`, with anything in the middle (the `*` character), as long as they also end in `fna.gz`, giving us:

```{bash, eval=TRUE}
ls ~/data/ref/bt2_indices/*fna.gz
```



command | output
-- | -----------
`sed` | takes text, from a file, or `piped` using the `|` command, and changes it
`-e` | `sed` has many options; we're using `-e` to tell it to use an `--expression` - i.e. regular expressions a.k.a. `regex`
`s/...` | we want `sed` to `s`ubstitute something for something different
`/.a./.b./` | `sed` will look between the first pair of slashes, and `s`ubstitute, swapping `a` for `b`
 `*` | allows any number of characters to be matched times, 0-infinity 
 `.` | will match any character - often put with `*` to give `.*`, meaning "match _any_ character, as often as possible" (see "[greedy]("https://en.wikipedia.org/wiki/Greedy_algorithm"))  
 `(...)` | text matched within the brackets in `/.a./` will be copied, and can then be pasted into `/.b./`. Multiple different things can be captured from `/.a./`, and pasted back into `/.b./` using `\1` for the first thing saved, `\2` for the second thing, `\3` for the third thing, etc.
 `\` | backslash 'escapes' the next character - `sed` will not try _match_ the next character, instead it interprets the next character's `regex` meaning - see [here](https://en.wikipedia.org/wiki/Regular_expression#POSIX_basic_and_extended) for a list.
`.../g` | here, we ask `sed` to do all of this `g`lobally, i.e. whereever it appears

So, the command `sed -e 's/.*\/\(.*\).fna.gz/\1/g'`:

command | output
-- | -----------
`sed` | calls `sed`
`-e` | activates regular expression mode
`s` | asks for substitutions
`/.*\/\(.*\)...` | part `a`: match as much as you can, i.e. `.*`, until you find the last `/` character; then save everything _after_ the last `/`, using the `( )`
`... .fna.gz` |  part `a`: this part is outside of the `( )` brackets above, so don't include it in the saved match from part `a`
`/\1/` | part `b`: paste in the bit saved from part `a`, i.e. everything after the last bracket

Is changed to the below, giving us a list of our filenames, without the `fna.gz` extension - really useful for calling inside scripts!

```{bash, eval=TRUE}
ls ~/data/ref/bt2_indices/*fna.gz | sed -e 's/.*\/\(.*\).fna.gz/\1/g'
```

---


<!-- To do this, we'll use another `for-loop`, which loops through the reference genomes we have, and builds an index for each one. It might be easier (to read, to look at, to type, etc.) if you rename your genomes to something simple, like `chicken.fna.gz` etc., but we won't take that step here. -->

<!-- # cow.853MB.fna.gz; 28m, 15GB rRAM, 9 threads, -->
<!-- time bowtie2-build --t 9 \ -->
<!--   ~/data/ref/bt2_indices/GCA_021234555.1_ARS-LIC_NZ_Jersey_genomic.fna.gz \ -->
<!--   ~/data/ref/bt2_indices/jersey_ARS-LIC.bt2 > \ -->
<!--   ~/data/ref/bt2_indices/jersey_ARS-LIC.bt2.buildlog -->

