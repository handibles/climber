##   kaiju   ======================

# https://kaiju.binf.ku.dk
# https://github.com/bioinformatics-centre/kaiju     # more likley		
# https://teaching.healthtech.dtu.dk/22126/index.php/Kaiju_exercise

## copy the DB to the existing /dev/shm/ (or mount & cp )
## alt : mount a memory device, and point DB to that : https://www.howtoforge.com/storing-files-directories-in-memory-with-tmpfs
# ls $KaDB/*{nodes,names,fmi}*
# TMP=/dev/shm
# time cp $KaDB/*{nodes,names,fmi}* $TMP/ & 

## likely initial DB mapping step also, substantial if nr_euk
## who knows what defaults actually mean:
#-e is allowed mismatches 	- def = 3
#-m is match length 		- def = 11 
#-s is match score 			- def = 65

## it was the tabs apparently
# KaDB=/data/databases/kaiju_v1.7.4/nr_euk/nr_euk     # lets try nr-euk! 

KaOUT=$WRK/4__kaiju
KaDB=/data/databases/kaiju/nr     # or not, its very big... 
mkdir $WRK/4__kaiju

module load kaiju

## test first  - 3m.
time kaiju \
    -v \
    -z 5 \
    -e 3 \
    -m 11 \
    -s 75 \
    -t $KaDB/nodes.dmp \
    -f $KaDB/kaiju_db_nr.fmi \
    -i $KDOUT/${TEST}_bt2decon_R1.fastq.gz \
    -j $KDOUT/${TEST}_bt2decon_R2.fastq.gz \
    -o $KaOUT/${TEST}_kaiju.out

## view - similar to bowtie2, but shows matches in amino acid seqs
less -S $KaOUT/${TEST}_kaiju.out

# what fraction assigned?
KA_C=$(grep -c -E '^C' $KaOUT/${TEST}*out) 
KA_U=$(grep -c -E '^U' $KaOUT/${TEST}*out) 

echo $(( $KA_C / $KA_U))


## same idea as kraken2:

# ```{bash, eval=FALSE}
echo '#!/bin/bash

#SBATCH --job-name=kaij_filt
#SBATCH --output=kaij_filt.%j.txt

#SBATCH --time=15:00

SAMPLE=$1
IN=$2
OUT=$3
DB=$4
THREADS=$5

time kaiju \
    -v \
    -z $THREADS \
    -e 3 \
    -m 11 \
    -s 65 \
    -t $DB/nodes.dmp \
    -f $DB/kaiju_db_nr.fmi \
    -i $IN/${TEST}_bt2decon_R1.fastq.gz \
    -j $IN/${TEST}_bt2decon_R2.fastq.gz \
    -o $OUT/${TEST}_kaiju__paired.out	


' > $MAT/slurm_kaij.sh  # SAMPLE IN OUT DB THREADS

```


```{bash, eval=FALSE}
# check first using $TEST
sbatch $MAT/slurm_kaij.sh $TEST $KDOUT $KROUT $K2REFDIR $KR_threads 

# pass all to slurm
for i in $( cat $MAT/samples );

  do sbatch $MAT/slurm_kaij.sh $i $KDOUT $KROUT $K2REFDIR $KR_threads ;

done
 
```

## taxonomy   ==================================================


```{bash, eval=FALSE}
##  issue solved by adding -l arg in kaiju2Table. However, (probably) still going to drop a lot of reads at the species level, so possibly necessary to compute a table for each depth 

ls -d $KaOUT/*/ | sed -E 's/.*\/(.*)\//\1/g' | sort | uniq | parallel -j 22 kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -r genus -l domain,phylum,class,order,family,genus -o $KaOUT/{}/{}_kaiju__summary_genus.tsv $KaOUT/{}/{}_kaiju__out.total
cat $KaOUT/*/*summary_genus.tsv | sed -E 's/.*\/(.*_kaiju__out.total.*$)/\1/' | sed -E 's/_kaiju__out.total//g' > $KaOUT/ms__paedcf_kaiju_genus.master

ls -d $KaOUT/*/ | sed -E 's/.*\/(.*)\//\1/g' | sort | uniq | parallel -j 22 kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -r phylum -l domain,phylum,class,order,family,phylum -o $KaOUT/{}/{}_kaiju__summary_phylum.tsv $KaOUT/{}/{}_kaiju__out.total
cat $KaOUT/*/*summary_phylum.tsv | sed -E 's/.*\/(.*_kaiju__out.total.*$)/\1/' | sed -E 's/_kaiju__out.total//g' > $KaOUT/ms__paedcf_kaiju_phylum.master

```


 # ===========================================================================================================

## tried stating memory needs (--mem) and changing to nr; tried stating ntask=1, ncpus = 5

# nano $MAT/slurms/slurm_kaiju.sh
# SAMPLE=$TEST
# IN=$KDOUT
# REFDIR=$TMP
# OUT=$KaOUT

## removed sbatch arg for --mem, specific nodes


# echo '
#!/bin/bash

# SBATCH –-job-name=kaiju_PUP
# SBATCH –-output=kaij_PUP
# SBATCH –-time=25:30
# SBATCH –-ntasks=1
# SBATCH –-cpus-per-task=5
# SBATCH –-mem=85gb
# SBATCH –-ntasks-per-node=1

## state yr var
SAMPLE=$1		
IN=$2
OUT=$3
REFDIR=$4


## make a kaiju box
mkdir $OUT/$SAMPLE
OUTDIR=$OUT/$SAMPLE


# ## kaiju on paired, sorted reads, avoid OOOError : seqtk sort  : https://bioinf.shenwei.me/seqkit/usage/#sort
# 
# # wayyyy too much e coli (databse representation) - increase score to tune the specificity of Greedy mode (no increase in mismsatches, -e)
# seqkit sort -n $IN/${SAMPLE}/${SAMPLE}_R1_kneaddata_paired_1.fastq.gz -o $OUTDIR/${SAMPLE}_R1_kneaddata_paired_1_stksort.fastq.gz
# seqkit sort -n $IN/${SAMPLE}/${SAMPLE}_R1_kneaddata_paired_2.fastq.gz -o $OUTDIR/${SAMPLE}_R1_kneaddata_paired_2_stksort.fastq.gz

time kaiju \
-v \
-z 5 \
-e 3 \
-m 11 \
-s 75 \
-t $REFDIR/nodes.dmp \
-f $REFDIR/kaiju_db_nr.fmi \
-i $OUTDIR/${SAMPLE}_R1_kneaddata_paired_1_stksort.fastq.gz \
-j $OUTDIR/${SAMPLE}_R1_kneaddata_paired_2_stksort.fastq.gz \
-o $OUTDIR/${SAMPLE}_kaiju__paired.out	>> $OUTDIR/${SAMPLE}.stout


## combine paired / unmatched
cat $OUTDIR/${SAMPLE}_kaiju__*.out > $OUTDIR/${SAMPLE}_kaiju__out.total

# formats exports etc  -  shoulda used species!  -  but catch tax here using -l 
if test -s $OUTDIR/${SAMPLE}_kaiju__out.total; then
kaiju2table -t $REFDIR/nodes.dmp -n $REFDIR/names.dmp -r species -l domain,phylum,class,order,family,genus,species -o $OUTDIR/${SAMPLE}_kaiju__summary.tsv $OUTDIR/${SAMPLE}_kaiju__out.total
else
  > $OUTDIR/${SAMPLE}.empty
fi


# ' > $MAT/slurms/slurm_kaiju.sh


## learm to slurm	
# $TEST
# ls $FILT/*R[12].fastq.gz | sed -e 's/.*\/\(.*\)_.*R..*/\1/g' | sort | uniq | grep SC143C1-TS_S35 | parallel -j 3 sbatch $MAT/slurms/slurm_kaiju.sh {} $KDOUT $KaOUT $KaDB
ls $FILT/*R[12].fastq.gz | sed -e 's/.*\/\(.*\)_.*R..*/\1/g' | sort | uniq | parallel -j 22 sbatch $MAT/slurms/slurm_kaiju.sh {} $KDOUT $KaOUT $KaDB


## catch pieces
# not sure why necessary. possibly can just re-fire all these missed parts at end. 
time kaiju \
-v \
-z 5 \
-e 3 \
-m 11 \
-s 75 \
-t $KaDB/nodes.dmp \
-f $KaDB/kaiju_db_nr.fmi \
-i $KaOUT/${TEST}_R1_kneaddata_unmatched_12.fastq.gz \
-o $KaOUT/${TEST}_kaiju__unmatched.out		>> $KaOUT/${TEST}.stout
# -i $OUTDIR/${SAMPLE}_R1_kneaddata_paired_1_stksort.fastq.gz \
# -j $OUTDIR/${SAMPLE}_R1_kneaddata_paired_2_stksort.fastq.gz \
# -o $OUTDIR/${SAMPLE}_kaiju__paired.out	>> $OUTDIR/${SAMPLE}.stout


## catch those inexplicably broken samples
ls $KaOUT/*/*empty | sed -E 's/.*kaiju\/(.*)\/.*/\1/g' | parallel -j 22 sbatch $MAT/slurms/slurm_kaiju.sh {} $KDOUT $KaOUT $KaDB



## taxonomy   ==================================================

##  issue solved by adding -l arg in kaiju2Table. However, (probably) still going to drop a lot of reads at the species level, so possibly necessary to compute a table for each depth 

ls -d $KaOUT/*/ | sed -E 's/.*\/(.*)\//\1/g' | sort | uniq | parallel -j 22 kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -r genus -l domain,phylum,class,order,family,genus -o $KaOUT/{}/{}_kaiju__summary_genus.tsv $KaOUT/{}/{}_kaiju__out.total
cat $KaOUT/*/*summary_genus.tsv | sed -E 's/.*\/(.*_kaiju__out.total.*$)/\1/' | sed -E 's/_kaiju__out.total//g' > $KaOUT/ms__paedcf_kaiju_genus.master

ls -d $KaOUT/*/ | sed -E 's/.*\/(.*)\//\1/g' | sort | uniq | parallel -j 22 kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -r phylum -l domain,phylum,class,order,family,phylum -o $KaOUT/{}/{}_kaiju__summary_phylum.tsv $KaOUT/{}/{}_kaiju__out.total
cat $KaOUT/*/*summary_phylum.tsv | sed -E 's/.*\/(.*_kaiju__out.total.*$)/\1/' | sed -E 's/_kaiju__out.total//g' > $KaOUT/ms__paedcf_kaiju_phylum.master


#     ## should probably put this within SLURM - gets correct NA for blank ranks, but lost at summary step
#     # parallel -j 22 kaiju-addTaxonNames -t $KaDB/nodes.dmp -n $KaDB/names.dmp -i {} -o {}.tax -r superkingdom,phylum,class,order,family,genus,species ::: $KaOUT/[L,S]*/*total 
#   
#   ## needs -p  ; still doesn't regularise the taxonomic output. 
#   ## however, note also that this will drop a load of sequences if unclear where they belong  -  maybe don't go with species?
#     parallel -j 22 kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -p -r species -o {}.tsv {} ::: $KaOUT/*/*tax
#     cat $KaOUT/*/*tax.tsv | sed -E 's/.*\/(.*_kaiju__pup.total.*$)/\1/' | sed -E 's/_kaiju__pup.total.tax//g' > $KaOUT/ms__paedcf_kaiju.master.tax
# 
# 
#     cat $KaOUT/*/*tax.tsv | sed -E 's/.*\/(.*_kaiju__pup.total).*/\1/' | head


## APPROACH  ====================================================

## ridiculously, have to then copy if to a local unix box, then USB it to here... 
scp -i ~/.ssh/id_rsa_ucc $KaOUT/*master* jamie@dunsh.ucc.ie:/home/jamie/
  ssh -i ~/.ssh/id_rsa_ucc jamie@dunsh.ucc.ie "scp -P 55555 ~/*master* jfg@143.239.154.14:/home/jfg/Dropbox/teag_sync"


# then from dunshire to asp

