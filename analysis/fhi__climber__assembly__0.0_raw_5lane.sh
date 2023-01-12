

## =============   re-nur  ================================================================

# databases for different tools
DB=/data/databases                                            # HPC only
HGR38_BT=/data/databases/hostremoval/Homo_sapiens/Bowtie2     # HPC only
K2REFDIR=/data/databases/kraken2                              # HPC only

# our overall data structure
RAW=/data/Food/primary/R0937_chickmicro
LANE=$RAW/lanemerge
WRK=/data/Food/analysis/R0602_microsupport/jamie.fitzgerald/r0937
BT_DB=$INES/ref

# scripts, and slurms for queueing up jobs
MAT=$WRK/Materials

# outputs etc
QC=$WRK/1__qc
FQC=$QC/fastqc
FILT=$WRK/2__filt
KDOUT=$WRK/3__knead
KROUT=$WRK/4__krak2

#!# set a sample-name to test things on
TEST=_N1_Nuria

mkdir $BT_DB $RAW $WRK $FILT $QC $FQC $KDOUT $KROUT $LANE $MAT
mkdir $FQC/fhi_michick_lane $FQC/fhi_michick_filt $FQC/fhi_michick_lane_multi $FQC/fhi_michick_filt_multi


## merge lanes1-4 and lane1  ================================================================

head $(find $RAW -name ${TEST}*R1*fastq.gz) #> $LANE/${TEST}_R1.txt
zcat $(find $RAW -name ${TEST}*R1*fastq.gz) > $LANE/${TEST}_R1.txt

# slow for days. with a subshell and not expansion - not sure why 
# par? nope, gzip stdin issues
$( find /data/Food/primary/R0937_chickmicro/ -name ${TEST}*R1*fastq.gz )
# R1
for i in $( cat $NUR/Materials/samples | sed 's/_S.*//' ); do zcat $( find /data/Food/primary/R0937_chickmicro/ -name ${i}*R1*fastq.gz ) > $LANE/${i}_merge_R1.fastq ; done
# R2
for i in $( cat $NUR/Materials/samples | sed 's/_S.*//' ); do zcat $( find /data/Food/primary/R0937_chickmicro/ -name ${i}*R2*fastq.gz ) > $LANE/${i}_merge_R2.fastq ; done
# compress
parallel -j 8 gzip {} ::: $RAW/lanemerge/*fastq


## FQC on the merged reads  ================================================================

module load multiqc

# write a slurm script first
echo '#!/bin/bash
#SBATCH --job-name=merge_fqc
#SBATCH --output=merge_fqc.%j.txt
#SBATCH --ntasks=35
#SBATCH --time=70:00

IN=$1
OUT=$2

module load fastqc multiqc
time fastqc -t 35 $IN/*fastq.gz -o $OUT
multiqc $OUT -o ${OUT}_multi
module unload fastqc multiqc
' > $MAT/slurm_fqc.sh

sbatch $MAT/slurm_fqc.sh $LANE $FQC/fhi_michick_lane
#multiqc $FQC/fhi_michick_lane -o $FQC/fhi_michick_lane_multi


## trimm qual and adapters  ================================================================


echo '>adapter_seq
CTGTCTCTTATACACATCT
>adapter_mate_seq
AGATGTGTATAAGAGACAG
>Transposase_Adap__for_tagmentation_1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Transposase_Adap__for_tagmentation_2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>PCR_primer_index_1
CAAGCAGAAGACGGCATACGAGATNNNNNNNGTCTCGTGGGCTCGG
>PCR_primer_index_2
AATGATACGGCGACCACCGAGATCTACACNNNNNTCGTCGGCAGCGTC
>PCR_primer_index_2_rc
GACGCTGCCGACGANNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT
>PCR_primer_index_1_rc
CCGAGCCCACGAGACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
>Transposase_Adap__for_tagmentation_2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>Transposase_Adap__for_tagmentation_1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>adapter_mate_seq_rc
CTGTCTCTTATACACATCT
>adapter_seq_rc
AGATGTGTATAAGAGACAG' > $MAT/fqc_trimmo_ill_ref.fa

# combine all those different parts!
ls $LANE/*fastq.gz | sed -e 's/.*\/\(.*\)_merge_R.*/\1/g' | sort | uniq > $MAT/samples


echo '#!/bin/bash

#SBATCH --job-name=bitcoinX
#SBATCH --output=trimmoRaw.txt
#SBATCH --open-mode=append
#SBATCH --ntasks=6
#SBATCH --time=32:00

SAMPLE=$1
IN=$2
OUT=$3
REFDIR=$4

echo ""
echo "   ## trimmomatic on ${SAMPLE} - CROP:145 HCROP:20 ILLWIND:2.30.10.5 SLID:6:15 MINL:100  --------------------  ## "
echo ""

module load trimmomatic
# trimmomatic - use backslash to separate rows
time trimmomatic PE \
$IN/${SAMPLE}_merge_R1.fastq.gz \
$IN/${SAMPLE}_merge_R2.fastq.gz \
$OUT/${SAMPLE}_R1_trimm.fastq.gz \
$OUT/${SAMPLE}_R1_trimm_unpaired.fastq.gz \
$OUT/${SAMPLE}_R2_trimm.fastq.gz \
$OUT/${SAMPLE}_R2_trimm_unpaired.fastq.gz \
CROP:145 \
HEADCROP:20 \
ILLUMINACLIP:$REFDIR/fqc_trimmo_ill_ref.fa:2:30:10:5 \
SLIDINGWINDOW:6:15 \
MINLEN:100 \
-threads 6 > $OUT/trimmo_${SAMPLE}.out 2>&1    # effective stout/err redirect

if [ ! -s $OUT/${SAMPLE}_R1_trimm.fastq.gz ] && [ ! -s $OUT/${SAMPLE}_R2_trimm.fastq.gz ]
then
	echo "   ## trimmo on ${SAMPLE} failed    ------------------    < ! > "
fi

module unload trimmomatic

' > $MAT/slurm_trimmo.sh

for i in $( cat $MAT/samples );
  do sbatch $MAT/slurm_trimmo.sh $i $LANE $FILT $MAT;
done

 
# check fqc to it
mkdir ${FILT}_unpaired
mv $FILT/*unpaired*fastq.gz ${FILT}_unpaired/
sbatch $MAT/slurm_fqc.sh $FILT $FQC/fhi_michick_filt

 
## make deco database  ================================================================
	#
	## already done, no need to redo
	#
	#module load bowtie2 samtools/1.10 parallel fastqc multiqc
	#
	## download your own genome (more or less) - 1.2 GB, 12min :: from the NCBI browser https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000001405.40/
	#curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/GCF_000001405.40/download?filename=GCF_000001405.40.zip" -H "Accept: application/zip"
	#
	## decoy genome - 9 MB
	#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5ss.fa.gz
	#
	#chmod u+x $MAT/trimmo_genomes_dl.sh
	#mkdir $RAW/ref
	#$MAT/trimmo_genomes_dl.sh $RAW/ref
	#
	#parallel gzip {} ::: $RAW/ref/*fna
	#
	#
	## number of processes/threads to use for jobs
	#BT_threads=10
	#BT_DB=/data/Food/analysis/R0936_michickeese/ref
	#mkdir $BT_DB
	#
	## build our BT2 database in slurm:
	#echo '#!/bin/bash
	#
	##SBATCH --job-name=uploadbitcoins
	##SBATCH --output=bt2DB.txt
	##SBATCH --ntasks=10
	##SBATCH --time=179:00
	#
	#INOUT=$1
	#
	## do check the names before you run! 
	#time bowtie2-build --large-index --threads 10 \
	#$INOUT/hs37d5ss.fa.gz \
	#$INOUT/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
	#$INOUT/deco.bt2
	#
	#' > $MAT/slurm_bt2db.sh
	#$MAT/slurm_bt2db.sh $BT_DB
	#
	#
## decontaminate FASTQs  ================================================================

echo '#!/bin/bash

#SBATCH --job-name=btcoinDL.%j
#SBATCH --output=deco_bt2.%j.txt
#SBATCH --ntasks=10

SAMPLE=$1
IN=$2
OUT=$3
REF=$4
THREADS=$5

# load necess
module load bowtie2 samtools/1.10

## stamp for the output
DAT_TIM=$(date +"%T %F")
echo "--------------------------------------------"
echo $DAT_TIM
echo sampl: $SAMPLE
echo input: $IN
echo outpt: $OUT
echo "--------------------------------------------"


echo "--   align   -------------------------------------"
# -x $REF/mult_ref \
time bowtie2 -p $THREADS \
-x $REF/deco.bt2 \
-1 $IN/${SAMPLE}_R1_trimm.fastq.gz \
-2 $IN/${SAMPLE}_R2_trimm.fastq.gz \
-S $OUT/${SAMPLE}_bt2_refmapped.sam


if [ $? = 0 ] && [ -s $OUT/${SAMPLE}_bt2_refmapped.sam ]
then
        echo "--   change format   -------------------------------------"
        time samtools view -bS $OUT/${SAMPLE}_bt2_refmapped.sam > $OUT/${SAMPLE}_bt2_refmapped.bam

        echo "--   remove matches   ------------------------------------"
        time samtools view -b -f 12 -F 256 $OUT/${SAMPLE}_bt2_refmapped.bam > $OUT/${SAMPLE}_bt2_decontamd.bam

		## we skip the sort step, as this aint genome assembly

        echo "--   to fastq, dump other   ------------------------------"
        time samtools fastq -@ $THREADS $OUT/${SAMPLE}_bt2_decontamd.bam -1 $OUT/${SAMPLE}_bt2decon_R1.fastq.gz -2 $OUT/${SAMPLE}_bt2decon_R2.fastq.gz -0 /dev/null -s /dev/null -n

        echo "---   $SAMPLE completed  :)   ----------------------------------"


else
        echo ""
        echo "< ! >---   align failed   ::   $SAMPLE   ---< ! >"
        echo ""
        echo ""

fi


## shrink intermediates
if [ -s $OUT/${SAMPLE}_bt2decon_R1.fastq.gz ] && [ -s $OUT/${SAMPLE}_bt2decon_R2.fastq.gz ] ;
then 
        echo "--   shrinking *AMs   ----------------------------------" ;
        gzip $OUT/${SAMPLE}_*.{bam,sam} &
fi

' > $MAT/slurm_deco_full.sh # SAMPLE  IN  OUT REF THREADS


# send just one sample to the script with "head -1"
for i in $( cat $MAT/samples );
    do sbatch $MAT/slurm_deco_full.sh $i $FILT $KDOUT $INES/ref 10 ;
done

## fix - not properly converted...
	#module load samtools/1.10
	#for i in $( cat $MAT/samples ); do samtools fastq -@ 10 $KDOUT/${i}_bt2_decontamd.bam -1 $KDOUT/${i}_bt2decon_R1.fastq.gz -2 $KDOUT/${i}_bt2decon_R2.fastq.gz -0 /dev/null -s /dev/null -n ; done
	#cat $MAT/samples | parallel -j 12 samtools fastq -@ 10 $KDOUT/{}_bt2_decontamd.bam -1 $KDOUT/{}_bt2decon_R1.fastq.gz -2 $KDOUT/{}_bt2decon_R2.fastq.gz -0 /dev/null -s /dev/null -n
	#parallel gzip {} ::: KDOUT/*am

# check outputs, check error messages


## kraken microbial Community ================================================================


# test, but also make for taxonomy
module load kraken2/2.1.1
KR_threads=5
time kraken2 --db $K2REFDIR \
$KDOUT/${TEST}_bt2decon_R1.fastq.gz \
$KDOUT/${TEST}_bt2decon_R2.fastq.gz \
--paired \
--confidence 0.15 \
--minimum-hit-groups 3 \
--minimum-base-quality 10 \
--report-zero-counts \
--gzip-compressed \
--threads $KR_threads \
--report $KROUT/${TEST}_test_kraken_report #\
#--unclassified-out $KROUT/${TEST}_test_kraken_unclass# \
#--output $KROUT/${TEST}_test_kraken_output



echo '#!/bin/bash

#SBATCH --job-name=bitcoinDB
#SBATCH --output=krak2_loop.%j.txt
#SBATCH --ntasks=10
#SBATCH --time=180:00

LIST=$1
IN=$2
OUT=$3
DB=$4
THREADS=$5

module load kraken2/2.1.1


# explicit about our parameters
CONF=0.15
MINH=3
MINQ=10


#one big loop / aw yeh aw yehh
for i in $( cat $LIST );
do
  if [ -s $OUT/${i}*_kraken_report ]
  then
    echo "## krak2 output for $i already in $OUT - stopping"
  else
    echo "## firing krak2 on $i ; conf $CONF: min-hits $MINH; min qual: $MINQ "

    time kraken2 --db $DB \
    $IN/${i}_bt2decon_R1.fastq.gz \
    $IN/${i}_bt2decon_R2.fastq.gz \
    --paired \
    --confidence $CONF \
    --minimum-hit-groups $MINH \
    --minimum-base-quality $MINQ \
    --threads $THREADS \
    --gzip-compressed \
    --report $OUT/${i}_kraken_report \
    --unclassified-out /dev/null/${i}_kraken_unclass# \
    --output /dev/null/${i}_kraken_output

  fi

done' > $MAT/slurm_krak2_loop.sh

# pass all to slurm
sbatch $MAT/slurm_krak2_loop.sh $MAT/samples $KDOUT $KROUT $K2REFDIR $KR_threads


## bracken estimates  ================================================================

### build data base first
#BR_kmer=35    # this is the default kmer length of the Kraken2 DB on the HPC
#BR_leng=125   # length post-trimm
#
### make a local install of bracken via github for building DB
#mkdir ~/bin ; cd ~/bin &&
#wget https://github.com/jenniferlu717/Bracken/archive/refs/heads/master.zip ; unzip master.zip && rm master.zip
#cd ~/bin/Bracken-master
#chmod +x install_bracken.sh ; ./install_bracken.sh
#
#
#echo '#!/bin/bash
#
##SBATCH --job-name=brack_db
##SBATCH --njobs=35
##SBATCH --output=brack_db.%j.txt
#
##SBATCH --time=240:00
#
#BRDB_DIR=$1
#BRDB_KMER=$2
#BRDB_LENG=$3
#BRDB_THREADS=$4
#
### not tested
#time ~/bin/Braken-master/bracken-build -d $BRDB_DIR -k $BRDB_KMER -l $BRDB_LENG -t $BRDB_THREADS
#
#' > $MAT/slurm_brak_db.sh # dir kmer length threads

## run on the existing, default kraken database
#sbatch slurm_brak_db.sh /data/databases/kraken2 $BR_kmer $BR_leng 20

module load braken
BR_r=125
BR_l=S
BR_t=10   # counts! not threads

for i in $(cat $MAT/samples );
  do bracken -d $K2REFDIR/ -i $KROUT/${i}_kraken_report -o $KROUT/${i}.bracken -r $BR_r -l $BR_l -t $BR_t ;
done                                                                                        

combine_bracken_outputs.py --files $KROUT/*.bracken -o $KROUT/krakenStnd_abundances.tsv

## install if necess
#mkdir $HOME/bin ; cd $HOME/bin ; git clone https://github.com/jenniferlu717/KrakenTools.git ; cd $HOME
~/bin/KrakenTools/kreport2mpa.py -r $KROUT/${TEST}_test_kraken_report -o $KROUT/${TEST}_kraken_mpa
grep -h '|s_' $KROUT/${TEST}_kraken_mpa | cut -f 1 | sort | uniq | sed 's/|/\t/g' > $KROUT/krakenStnd_taxonomy.tsv

less -S $KROUT/krakenStnd_taxonomy.tsv
less -S $KROUT/krakenStnd_abundances.tsv

## also, get the combined kraken reports, for procrustes of community
~/bin/KrakenTools/combine_kreports.py -r $KROUT/*Nuria_kraken_report -o $KROUT/krakenStnd_kreportCombo.tsv

less -S $KROUT/krakenStnd_kreportCombo.tsv
# table of sorts
grep -E 'perc\stot_all' $KROUT/krakenStnd_kreportCombo.tsv | sed 's/#//' > $KROUT/krakenStnd_kreportCombo_count.tsv
grep -E '\sS1?\s' $KROUT/krakenStnd_kreportCombo.tsv >> $KROUT/krakenStnd_kreportCombo_count.tsv

## in R  -  



## Kaiju microbial Community ================================================================

KaDB=/data/databases/kaiju/nr     # nr-euk very big..
Ka_THREADS=10
KaOUT=$WRK/4__kaiju
mkdir $KaOUT
# ulimit -s unlimited    # nasa mem trick


## slurm-loop  kaiju   =======================================

echo '#!/bin/bash

#SBATCH --job-name=NFTcompil
#SBATCH --output=kaij_loop.%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=150:00

LIST=$1
IN=$2
OUT=$3
DB=$4
THREADS=$5


## one big slurm to for-loop over kaiju, thereby avoiding mem issues

# operate over some list of sample names
for i in $( cat $LIST );
do
  if [ -s $OUT/${i}*_kaiku ]
  then
    echo "## kaiju output for $i already in $OUT - stopping"
  else
    echo "## firing kaiju on $i"

    module load kaiju/1.7.4
    time kaiju \
    -v \
    -e 3 \
    -m 11 \
    -s 65 \
    -z $THREADS \
    -t $DB/nodes.dmp \
    -f $DB/kaiju_db_nr.fmi \
    -i $IN/${i}_bt2decon_R1.fastq.gz \
    -j $IN/${i}_bt2decon_R2.fastq.gz \
    -o $OUT/${i}_kaiju

    # count things
    if [ -s $OUT/${i}*_kaiju ]
    then
      KAI_C=$(grep -cE '^C' $OUT/${i}_kaiju ) &
      KAI_U=$(grep -cE '^U' $OUT/${i}_kaiju ) 
      KAI_TOT=$(echo "$KAI_C + $KAI_U" | bc)
      KAI_PC=$(echo "scale=2; ($KAI_C / $KAI_TOT)*100" | bc)
      echo "## kaiju sample processed: ${i} : total classified: ${KAI_PC}% (total: $KAI_TOT read-pairs)  ------"
    else
      echo "## no output - kaiju for sample ${i} failed"
    fi

  fi
    
done


' > $MAT/slurm_kaij_loop.sh  # **LIST** IN OUT DB THREADS


## fire all in serial
sbatch $MAT/slurm_kaij_loop.sh  $MAT/samples $KDOUT $KaOUT $KaDB $Ka_THREADS ;


for i in $( cat $MAT/samples );
do 
  KAI_C=$(grep -cE '^C' $KaOUT/${i}_kaiju )
  KAI_U=$(grep -cE '^U' $KaOUT/${i}_kaiju ) 
  KAI_TOT=$(echo "$KAI_C + $KAI_U" | bc)
  KAI_PC=$(echo "scale=2; ($KAI_C / $KAI_TOT)*100" | bc)
  echo "## kaiju sample processed: ${i} : total classified: ${KAI_PC}% (total: $KAI_TOT read pairs)  ------"
done


## all taxonomy together, *AT THE SPECIES LEVEL* - some issue solved by adding -l arg in kaiju2Table
kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -r species -c 10 -l domain,phylum,class,order,family,genus,species -o $KaOUT/kaiju_summary_species.tsv $KaOUT/*_kaiju

