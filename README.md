# `CLI`, `M`icro`b`ial `E`cology, and `R`

Some basic steps in microbial ecology, focusing on the processing of `2ndGen` Illumina `fastq` data, into either `amplicon` (e.g. 16S) or `metagenomic` (e.g. shotgun) datasets, followed by ecology-based analysis of the communities and patterns we find in that data.

```mermaid
  graph TD;
      A-->B;
      A-->C;
      B-->D;
      C-->D;
```



```mermaid

flowchart LR
    fastq1[(sequencing run\ndata)] ==> |download data,\npresumably from illumina| qc1(((fa:fa-arrow-up-right-from-square 1-3: clean\nsequences)));
    qc1 ==> tax1(((fa:fa-arrow-up-right-from-square 4: taxonomy\nassignment)));
    qc1 ==> func1(((fa:fa-link-slash X: Function)));
    qc1 -.-> mag1(((fa:fa-link-slash Y: Metagenomic\nassembly)));
    mag1 -.- |figuring this out!| mag1;
    tax1 ==> comm1;
    func1 ==> comm1;
    comm1 -.- |fa:fa-ban avoid endless\nrecursion\nif at all\npossible...| comm1;
    comm1(((fa:fa-arrow-up-right-from-square microbial\necology)));
    comm1 --> R2(diversity);
    comm1 --> R3(abundance);
    comm1 --> R4(association);
    R2 --> R5(visualise);
    R3 --> R5;
    R4 --> R5;
    click qc1 "https://handibles.github.io/climber/documents/shotgun_assembly.html#2_-_quality_control_(_check)";
    click tax1 "https://handibles.github.io/climber/documents/shotgun_assembly.html#4_-_Microbiome_Community_Profiling";
    click comm1 "https://handibles.github.io/climber/documents/data_to_R.html";
    click func1 "https://teagasc.ie/unknowable_truths";
    click mag1 "https://teagasc.ie/unknowable_truths";

```


### Metagenomic data (i.e. shotgun)

  * <a href="https://handibles.github.io/climber/documents/shotgun_assembly.html">`-->` Jump straight to the tutorial on metagenomic shotgun assembly </a> 

As above, the tutorial covers the following steps:

  0. Setting up your analysis - `bash` and friends
  1. Checking your sequence data - `FastQC` & `MultiQC`
  2. Sequencing QC - filtering and trimming your sequences - `Trimmomatic`
  3. Sequencing QC - purifying your sequences - `BowTie2`
  4. Metagenomic Community profiling - `Kraken2` & `Bracken`


We also move through <a href="documents/data_to_R.html">importing output from `Kaiju` or `Kraken2+Bracken` into `R`: </a>.

  * importing data into `R` - generating a count matrix, taxonomic table, and phyloseq object from metagenomic data


This metagenomic workflow is also present in simple, no-nonsense, `raw code` (note there might be differences to the complete workflow above).

  * <a href="documents/shotgun_assembly_raw.html">`raw code only of the metagenomic shotgun assembly`</a> - as above, less explanation


### Amplicon data (e.g. 16S)

Forthcoming. The initial steps (setup, get data, QC) are very similar in most cases (remember to cut off your primers!), but are followed by a denoising step (`DADA2`) and optionally an attempt to predict the metabolic capabilities of the communities at hand (`PICRUSt2`).


### Microbial Ecology (and `R`)

This is the real magic, and we get to make _pictures_. This part might be updated as annotated code ahead of actual tutorials - it's a massive collection of huge topics..


<a href="documents/mb6302__preamble.html">`see also here`</a>

---

This guide to metagenomic analysis continues to be updated (April, 3023). All (+/-)feedback is welcome: simply throw objects/comments directly at me, or [drop us a line at the related repo](https://github.com/handibles/climber/issues). <a href="documents/climber_todo.html">`</a>

all the best!  

<a href="https://www.fhi.ie/project/jamie-fitzgerald/"> `Jamie` </a>

<img src="vis/ucc.png" width="150" align="center" /> | <img src="vis/teag.png" width="150" align="center" /> | <img src="vis/v1group.png" width="150" align="center"/>
