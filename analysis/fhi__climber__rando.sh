

## make corresponding dirs

# databases for different tools
DB=/data/databases                                            # HPC only
HGR38_BT=/data/databases/hostremoval/Homo_sapiens/Bowtie2     # HPC only
K2REFDIR=/data/databases/Kraken2_GTDB_r89_54k                 # HPC only
BT_DB=/data/Food/analysis/R0936_redcheese/ref

# our overall data structure
RAW=/data/Food/primary/R0936_redcheese
INES=/data/Food/analysis/R0936_redcheese
WRK=/data/Food/analysis/R0602_microsupport/r0936

# scripts, and slurms for queueing up jobs
MAT=$WRK/Materials
SLURM=$WRK/slurms

# outputs etc
QC=$INES/1__qc
FQC=$QC/fastqc
FILT=$INES/2__filt
KDOUT=$WRK/3__knead
KROUT=$WRK/4__krak2

#!# set a filename to test things on
TEST=Q4T8

module load parallel fastqc multiqc trimmomatic bowtie2 kraken2 braken samtools/1.10

mkdir $WRK $KDOUT $KROUT $MAT


## skipping to the decontam stage


## deco scripts


# Bowtie2 align sequences - 13m on sample Q4T8 (1.1GB F, 1.1GB R fastq.gz)
echo '#!/bin/bash

#SBATCH --job-name=deco_bt2align
#SBATCH --output=decoalign.txt
#SBATCH --open-mode=append

SAMPLE=$1
IN=$2
OUT=$3
REF=$4
THREADS=$5

# -x $REF/mult_ref \
time bowtie2 -p $THREADS \
-x $REF/deco.bt2 \
-1 $IN/${SAMPLE}_R1_trimm.fastq.gz \
-2 $IN/${SAMPLE}_R2_trimm.fastq.gz \
-S $OUT/${SAMPLE}_bt2_refmapped.sam

' > $MAT/slurm_deco_bt2align.sh # SAMPLE  IN  OUT REF THREADS


# convert bowtie2 SAM output to BAM -  10m on sample Q4T8, 2GB output
echo '#!/bin/bash

#SBATCH --job-name=deco_sambam
#SBATCH --output=deco_sambam.txt
#SBATCH --open-mode=append

SAMPLE=$1
IN=$2
OUT=$3
REF=$4
THREADS=$5

## change format
time samtools view -bS $OUT/${SAMPLE}_bt2_refmapped.sam > $OUT/${SAMPLE}_bt2_refmapped.bam
# if [ -s $OUT/${SAMPLE}_bt2_refmapped.bam & ] ; rm $OUT/${SAMPLE}_bt2_refmapped.sam &

' > $MAT/slurm_deco_sambam.sh # SAMPLE  IN  OUT REF THREADS


# filter contaminants from BAM -  5m on sample Q4T8
echo '#!/bin/bash

#SBATCH --job-name=deco_bamfilt
#SBATCH --output=deco_bamfilt.txt
#SBATCH --open-mode=append

SAMPLE=$1
IN=$2
OUT=$3
REF=$4
THREADS=$5

# samtools removes reads which match (contaminant) reference genome
time samtools view -b -f 12 -F 256 $OUT/${SAMPLE}_bt2_refmapped.bam > $OUT/${SAMPLE}_bt2_decontamd.bam
#rm $OUT/${SAMPLE}_bt2_refmapped.bam &

' > $MAT/slurm_deco_bamfilt.sh # SAMPLE  IN  OUT REF THREADS


# resort BAM to sensible order
echo '#!/bin/bash

#SBATCH --job-name=deco_bamsort
#SBATCH --output=deco_bamsort.txt
#SBATCH --open-mode=append

SAMPLE=$1
IN=$2
OUT=$3
REF=$4
THREADS=$5

# samtools re-roganises reads 
time samtools sort -n -m 5G -@ $THREADS $OUT/${SAMPLE}_bt2_decontamd.bam -o $OUT/${SAMPLE}_bt2_decontamd_sort.bam
#rm $OUT/${SAMPLE}_bt2_decontamd.bam &

' > $MAT/slurm_deco_bamsort.sh # SAMPLE  IN  OUT REF THREADS


# convert decontaminated, sorted BAM, to FASTQ
echo '#!/bin/bash

#SBATCH --job-name=deco_bamfq
#SBATCH --output=deco_bamfq.txt
#SBATCH --open-mode=append

SAMPLE=$1
IN=$2
OUT=$3
REF=$4
THREADS=$5

# stream outputs of non-interest to /dev/null
time samtools fastq -@ $THREADS $OUT/${SAMPLE}_bt2_decontamd_sort.bam -1 $OUT/${SAMPLE}_bt2decon_R1.fastq.gz -2 $OUT/${SAMPLE}_bt2decon_R2.fastq.gz -0 /dev/null -s /dev/null -n

#if [ -f $OUT/${SAMPLE}_bt2decon_R2.fastq.gz & -f $OUT/${SAMPLE}_bt2decon_R2.fastq.gz ] ; then rm $OUT/${SAMPLE}_*.[bam,sam] & ; fi

' > $MAT/slurm_deco_bamfq.sh # SAMPLE  IN  OUT REF THREADS


# ultimate test
sbatch $MAT/slurm_deco_bt2align.sh $TEST $FILT $KDOUT $BT_DB 10   # 0.04%, 12min

sbatch $MAT/slurm_deco_sambam.sh $TEST $FILT $KDOUT $BT_DB 10     # "no" STOUT, 5m

sbatch $MAT/slurm_deco_bamfilt.sh $TEST $FILT $KDOUT $BT_DB 10    # "no" STOUT, 6m

# sbatch $MAT/slurm_deco_bamsort.sh $TEST $FILT $KDOUT $BT_DB 10    # 5m  - works
sbatch $INES/Materials/slurm_deco_bamsort.sh $TEST $FILT $KDOUT $BT_DB 10    # ...

sbatch $MAT/slurm_deco_bamfq.sh $TEST $FILT $KDOUT $BT_DB 10      # 2m



## full script

echo '#!/bin/bash

#SBATCH --job-name=deco_bt2.%j
#SBATCH --output=deco_bt2.%j.txt
#SBATCH --open-mode=append

SAMPLE=$1
IN=$2
OUT=$3
REF=$4
THREADS=$5


## stamp for the output
DAT_TIM=$(date +"%T %F")
echo "----------------------"
echo $DAT_TIM 
echo sampl: $SAMPLE 
echo input: $IN
echo outpt: $OUT
echo "----------------------"


# -x $REF/mult_ref \
time bowtie2 -p $THREADS \
-x $REF/deco.bt2 \
-1 $IN/${SAMPLE}_R1_trimm.fastq.gz \
-2 $IN/${SAMPLE}_R2_trimm.fastq.gz \
-S $OUT/${SAMPLE}_bt2_refmapped.sam


## change format
time samtools view -bS $OUT/${SAMPLE}_bt2_refmapped.sam > $OUT/${SAMPLE}_bt2_refmapped.bam


# samtools removes reads which match (contaminant) reference genome
time samtools view -b -f 12 -F 256 $OUT/${SAMPLE}_bt2_refmapped.bam > $OUT/${SAMPLE}_bt2_decontamd.bam


# samtools re-roganises reads 
time samtools sort -n -m 5G -@ $THREADS $OUT/${SAMPLE}_bt2_decontamd.bam -o $OUT/${SAMPLE}_bt2_decontamd_sort.bam


# stream outputs of non-interest to /dev/null
time samtools fastq -@ $THREADS $OUT/${SAMPLE}_bt2_decontamd_sort.bam -1 $OUT/${SAMPLE}_bt2decon_R1.fastq.gz -2 $OUT/${SAMPLE}_bt2decon_R2.fastq.gz -0 /dev/null -s /dev/null -n

## reinstate when tested properly
#if [ -s $OUT/${SAMPLE}_bt2decon_R1.fastq.gz & -s $OUT/${SAMPLE}_bt2decon_R2.fastq.gz ] ; then rm $OUT/${SAMPLE}_*.[bam,sam] & ; fi


' > $MAT/slurm_deco_full.sh # SAMPLE  IN  OUT REF THREADS

sbatch $MAT/slurm_deco_full.sh $TEST $FILT $KDOUT $BT_DB 10


## real problem shooting   ==============================================================
    # 
    # # make a new test
    # zcat $FILT/Q6T10_R1_trimm.fastq.gz | head -4000 > $KDOUT/test_R1_trimm.fastq
    # zcat $FILT/Q6T10_R2_trimm.fastq.gz | head -4000 > $KDOUT/test_R2_trimm.fastq
    # gzip $KDOUT/test*fastq
    # 
    # TOY=test
    # 
    # sbatch $MAT/slurm_deco_bt2align.sh $TOY $KDOUT $KDOUT $BT_DB 10   # 0.04%, 12min
    # sbatch $MAT/slurm_deco_sambam.sh $TOY $FILT $KDOUT $BT_DB 10     # no STOUT, 5m
    # sbatch $MAT/slurm_deco_bamfilt.sh $TOY $FILT $KDOUT $BT_DB 10    # ... ?
    # sbatch $MAT/slurm_deco_bamsort.sh $TOY $FILT $KDOUT $BT_DB 10    # 5m
    # sbatch $MAT/slurm_deco_bamfq.sh $TOY $FILT $KDOUT $BT_DB 10
    # 
    # 
    # sbatch $MAT/slurm_deco_full.sh $TEST $KDOUT $KDOUT $BT_DB 10   # 0.04%, 12min




