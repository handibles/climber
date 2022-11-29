
## metageomic assemblies    ========================================================== ##

  # metahit for metagenomic re-assembly    ***
  # metaSPAdes for metagenomic re-assembly...
  # unicycler for assembling from isolates
  # tricycler for consensus assembly of genomes

  ## via SM:
  # https://github.com/rrwick/Unicycler
  # https://github.com/rrwick/Trycycler/wiki/Guide-to-bacterial-genome-assembly) 
  # https://f1000research.com/articles/8-2138  F1000 long-read assembly comparison

## ----------------------------------------------------------------------------------- ##
## =================================================================================== ##


##   s e t u p   ## --------------------------------------------------------

REF=/data/databases                                            # HPC only
SEQ=~/basemount/Projects/HarshMicroMatrixExp23_2022-11-09T16_17_33_7f06fdd/AppSessions.v1/HarshMicroMatrixExp23/Properties/Output.BiologicalSamples
# our overall data structure
RAW=/data/Food/primary/R0602_microsupport/jamie.fitzgerald/fhi__mimat
WRK=/data/Food/analysis/R0602_microsupport/jamie.fitzgerald/fhi__mimat

# scripts, and slurms for queueing up jobs
MAT=$WRK/Materials

# outputs etc
QC=$WRK/1__qc
FILT=$WRK/2__filt
DECO=$WRK/3__knead
KRAK=$WRK/4__krak2
KAIJ=$WRK/4__kaiju

#!# set a sample-name to test things on
TEST=FHI9_S9


##   m e g a h i t   ## --------------------------------------------------------

echo '#!/bin/bash

#SBATCH --job-name=megahit.kstep10
#SBATCH --output=megahit.kstep10.%j.slrm
#SBATCH --ntasks=12
#SBATCH --time=300:00

SAMPLE=$1
IN=$2
OUT=$3

module load megahit
megahit -1  $IN/${SAMPLE}_bt2decon_R1.fastq.gz -2 $IN/${SAMPLE}_bt2decon_R2.fastq.gz -t 12 --k-step 10 -o ${OUT} --out-prefix ${SAMPLE}_

' > $MAT/slurm_megahit.kstep10.sh # SAMPLE IN OUT  -  feel free to use a better name too

GENO=$WRK/6__genomes
mkdir $GENO

# test
do sbatch $MAT/slurm_megahit.kstep10.sh $TEST $KDOUT $GENO/test

# adding stupid folders for stupid reasons
for i in $( ls $MAT/samples )
do sbatch $MAT/slurm_megahit.kstep10.sh $i $KDOUT $GENO/full
done



