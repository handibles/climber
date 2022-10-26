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


<!-- # --- -->
<!-- # title: 'Red Cheese - Preamble' -->
<!-- # author: "jfg / ic" -->
<!-- # date: "`r format(Sys.time(), '%d %B, %Y')`" -->
<!-- # output: -->
<!-- #   html_document: -->
<!-- #     toc: TRUE -->
<!-- #   pdf_document:  -->
<!-- #     toc: TRUE -->
<!-- # --- -->

See [RMD reference](https://rmarkdown.rstudio.com/docs/reference/html_document.html) for all header terms

  0. <a href="/climber/documents/0.setup.html">Setting up your analysis</a>
  1. <a href="/climber/documents/1.checkdata.html">Checking your sequence data</a>
  2. <a href="/climber/documents/2.seqcleaning.html">Sequencing QC - filtering and trimming your sequences</a>
  3. <a href="/climber/documents/3.seqpurity.html">Sequencing QC - purifying your sequences - **not finished**</a>

screen -list
sbatch 
squeue | less

mv fa fna


catching failures from logs, compare etc. 

```
# 3 - decontaminate data (& check)
@jfg: there is, in addition, a database of human sequences which are poorly represented in the human genome, but are highlighted in Breitweizer et al.'s 2019 paper ["Human contamination in bacterial genomes has created thousands of spurious proteins"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6581058/). Currently looking for these sequences but not found yet. possibly pull out of kneaddata?
```

no paramters are explored for bt build align, or samtools. Samtools format will require its own giude...

rename, and standardising the name of samples & extentions etc.

# to do
```{bash, eval = FALSE}
 # intro on bash tips
    #top, htop
    #ls, ls -lash
    #mkdir
    #cd
    #rm
    #rmdir
    #mv
    #echo, cat
    #intro on scp, sFTP, rsync
    #wc 
    #sed
    # time
    #hup
    #comm

 # intro on the HPC

 # intro on slurm

 # intro on R tips

 # intro on screen

 # intro on \

 # intro on |
 
 # intro on >
 
 # intro on globbing
 
 # intro on loops
 
 # intro on parallel
 
 # intro on shell exapansion
 
 # intro on rsync and backups

 # how to read an error

 # how to make this fake fastq data

 # reference material
 
 # bibliography

```

---

# Introduction

#### What this is, and who it is for 

This document is for someone (anyone!) trying to do the following:

  * ...
  * ...
  * ...

#### What this document is **not**

#### When something goes wrong:

### For **help** with a program:


---

# 0 setup your environment

## the HPC

The HPC (high performance compute) cluster is a system of 12 high powered computers (```nodes```) connected to a main access computer (```head```). We use ```ssh``` to access the HPC's head, and then **we always transfer to a node**. If you do not do this, it causes big problems, and you will be punished by magic (your user will be locked, any programs you run will be cancelled, and your project will be stuck). 

The HPC is just some really big computers that share one huge hard drive. For more, better info on Teagasc's HPC, see [this link](http://hcux400.teagasc.net) (note you must be connected to the teagasc HPC to access it). There is also [a lot of general info on the web](https://duckduckgo.com/?q=intro+hpc+slurm&t=vivaldi&atb=v314-1&ia=web), as many places have different HPCs. 


```{bash,eval=FALSE}
    # ssh/putty arrives here at the 'head'
    [[______________hcux400______________]]    
                      |
                      |   # we   A L W A Y S   change to a node!!
                      |
                      V
    # then we ssh AGAIN, to the nodes (1 & 11 are broken)
    [ computeXX ][ compute02 ][ compute03 ]
    
    [ compute04 ][ compute05 ][ compute06 ]
    
    [ compute07 ][ compute08 ][ compute09 ]
    
    [ compute10 ][ computeXX ][ compute12 ]

```


## connecting via SSH

To connect to computers remotely, we use a tool called ```ssh``` ([secure shell protocol](https://en.wikipedia.org/wiki/Secure_Shell)). While many computers (e.g. ubuntu) come with ```ssh``` software installed. Windows needs to add a program to provide ```ssh``` functions. To use the Teagasc HPC, you need to connect to the Teagasc network: if on site, this is done already. If remote/wifi, log in to the Teagasc VPN. Next:

 * **windows** : download [```putty``` (standalone version, x86)](www.putty.org/downloads) from the internet, and put it somewhere useful.
 * **windows** : click to open it - enter ```teagaschpc``` as hostname/IP address, and use port ```22```, which will be the default - then click ```Open```!
 * **mac/linux** computer : try opening a ```terminal``` window, and entering the same 

When you connect to the teagasc HPC with ```putty```, you are asked for:

  * your HPC ```username``` - Teagasc IT should have provided you with this
  * your HPC ```password``` - don't worry if you can't see it when you type!

```{bash, eval = FALSE}
# this is approximate! it might not look exactly like this
(base) your.name@hcux400:~$ 
```
  
The window should look like the above (system dependent), which is the start of the journey - you're now working at the ```bash``` commandline! Note the following pieces:

  * ```(base)``` - this part indicates an active ```conda environment```, which gives you extra powers, depending on what has been automatically "activated" when you login
  * ```your.name``` - this is you, the user
  * ```hcux400``` - this is the part of the hpc your using to run programs. ```hcux400``` is the head node.
  * ```~``` - this is the current working dir(ectory), where you are doing work and making files. "```~```" is the shortcut to a user's ```home``` dir, i.e. ```/home/your.name```.
  
_Note:_ Many other ```ssh``` methods are out there, but aren't discussed here:

 * [cygwin](https://www.cygwin.com/install.html) is like adding a mini-linux computer to your computer, and can be used for many different things:
 * [mobixterm](https://mobaxterm.mobatek.net/download.html) I have note used.
 * if you are on Linux and don't have ```ssh```, try ```apt-get install open-ssh``` (you'll need admin privileges).
 * **NB:** [RStudio](https://www.rstudio.com/products/rstudio/#rstudio-desktop) has a ```terminal``` - you may be able to run your ```ssh``` sessions out of RStudio, instead of using a ```terminal``` standalone or ```putty``` etc.

The first time we connect to the HPC with ```putty``` using ```ssh```, we set a safe, strong password, using ```passwd``` and following the instructions. The username and password will start the same as those used to log in to the HPC using putty/ssh. The next time, we will have to use the new password that we create below: 

```{bash, eval = FALSE}

login:
  user
  password
  
change password:
  passwd

```

  
  
### change node

```{bash, eval = FALSE}
ssh compute05
  
  ## set up shortcut

  which sftp
  which scp


## check the data partition

  /home/ines.calvete  = ~
  
  /data/Food/Primary/R0936_redcheese/
  
  /data/Food/Analysis/R0936_redcheese/
  

## check the analysis partition



## set your dir variables

## set your assumptions

## load your modules

```

---

## Remember:

Pretend that I am saying the following things very loudly at you, and giving you gestures of great importance:

  * you **must** be connected to the Teagasc network/VPN in order to access the  HPC.**
  * when you login to the HPC, you must **always** change to a ```node``` before you do any work 
  

## where does data go?

## basespace or similar

## some osrt of VPN:

The forticlient debacle

## ssh or fpt


---

# Reading / Reference

  * [amazing thread on the effects of trimming your data - from 900K to 600K kmers, as trimming removes fake kmer patterns](https://www.biostars.org/p/212136/). Go on Brian Bushnell!
  
  

---


# Second Document :: the handy commands

```{bash, eval = FALSE}
pwd

ls
ls .
ls /home/ines/calvete
ls ~

ls /data/Food/analysis/R0936_redcheese
ls /data/Food/primary/R0936_redcheese

cd 

nano


```


## declaring variables

```{bash, eval = FALSE}
ANA=/data/Food/analysis/R0936_redcheese
cd $ANA
```


## setting loops

## running parallel

## first functions

## get data

### qiita; qiita search terms; download

### ENA; search terms; download

	red cheese metagenomic illumina
	red cheese
	cheese # start broad!
   > projewcts :
	PRJNA295825	raw sequence reads from analysis of red-smear-cheese
The cheese smear of a Swiss semi-hard red-smear cheese variety was analyzed for comparison of the microbial surface smear community patterns of cheeses differing in their smear quality.


Here, we are taking a risk picking a random study - we could get anything back. However, the sequencing was carried out by GATC (a major company in the sequencing space), so should be relatively sane as long as the samples were extracted in a sensible way

In fact, when we click on the "Show Column Selection" tab, showing uns the available metadata, and enable instrument model, we see that this is all 454 pyrosequencing data (from the era before illumina). Nix'd!

### moving and un/zipping data



