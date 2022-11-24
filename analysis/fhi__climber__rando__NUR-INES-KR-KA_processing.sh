

  walkthrouhg for contaminated NURIA data

## todo
rearrange NUria samples
what database

INES: 
kai - ....
kr2 - done

NUR:
kai - done
kr2 - ...download, resort


## to show
comparison or evaluation of different db/methods
show what the assignment levels are 
	# R/Python is different,		
	# bash is different,			
	# HPC is a different environment

## done
two ines samples - processed
put samples in share	
 - ${INES}.${NUR} phylo basic Kaiju
 - ${INES}.${NUR} phylo basic Kraken2


## mondayt 21.11.22

# kraj2_stnd output in $INES/4__krak2
# redoing to be sure that we're using the rihgt db with krak2/bracken
mkdir $MSDAT/r0936/4__krak2_default
KR_threads=5
sbatch $MAT/slurm_krak2_loop.sh $INES/Materials/samples $MSDAT/r0936/3__knead $MSDAT/r0936/4__krak2_default $K2REFDIR $KR_threads

BR_r=130
BR_l=S
BR_t=10   # counts! not threads

module load braken

KROUT=$SHAR/r0936/4__krak2
for i in $(cat $INES/Materials/samples );
  do bracken -d $MSDAT/ref/ -i ${KROUT}_default/${i}_kraken_report -o ${KROUT}_default/${i}.bracken -r $BR_r -l $BR_l -t $BR_t ;
done                                                                                        

## combine reports to one monolith.tsv
combine_bracken_outputs.py --files ${KROUT}_default/*.bracken -o ${KROUT}_default/fhi_redch_krakenStnd_abundances.tsv

## remake the def test report, --report-zero-counts
	# as in kraken2 part

~/bin/KrakenTools/kreport2mpa.py -r ${KROUT}_default/${TEST}_test_kraken_report -o ${KROUT}_default/${TEST}_kraken_mpa
 
grep -h '|s_' ${KROUT}_default/${TEST}_kraken_mpa | cut -f 1 | sort | uniq  sed 's/|/\t/g' > ${KROUT}_default/fhi_redch_krakenStnd_taxonomy.tsv




## ====================
link for Kraken2  
  
graphical representations!!	


# INES  ==============================================

## refire Q4T7
KaDB=/data/databases/kaiju/nr
WRK=$MSDAT/r0936
KDOUT=$WRK/3__knead
KaOUT=$WRK/4__kaiju
KROUT=$MSDAT/r0936/4__krak2_default
# mkdir $KaOUT
TEST=Q4T7

# deco - skipped fastq day
module load samtools/1.10
gzip -d $KDOUT/${TEST}_bt2_decontamd.bam.gz
samtools fastq -@ 8 $KDOUT/${TEST}_bt2_decontamd.bam -1 $KDOUT/${TEST}_bt2decon_R1.fastq.gz -2 $KDOUT/${TEST}_bt2decon_R2.fastq.gz -0 /dev/null -s /dev/null -n
gzip $KDOUT/${TEST}_bt2_decontamd.bam &


# kraken2
module load kraken2/2.1.1
KR_threads=5
K2REFDIR=$MSDAT/ref   # 130mers

time kraken2 --db $K2REFDIR \
$KDOUT/${TEST}_bt2decon_R1.fastq.gz \
$KDOUT/${TEST}_bt2decon_R2.fastq.gz \
--paired \
--confidence 0.5 \
--minimum-hit-groups 5 \
--minimum-base-quality 20 \
--gzip-compressed \
--threads $KR_threads \
--report $KROUT/${TEST}_kraken_report \
--output $KROUT/${TEST}_kraken_output

module load braken
BR_r=130
BR_l=S
BR_t=10   # counts! not threads
bracken -d $K2REFDIR/ -i $KROUT/${TEST}_kraken_report -o $KROUT/${TEST}.bracken -r $BR_r -l $BR_l -t $BR_t
lk $KROUT/*.bracken
combine_bracken_outputs.py --files $KROUT/*.bracken -o $KROUT/krakenStnd_abundances.tsv


## redo kaiju for $TEST also
KaDB=/data/databases/kaiju/nr
ls -lsh $MSDAT/r0936/4__kaiju

# installed via conda 
conda activate kaijamie

KaOUT=$MSDAT/r0936/4__kaiju
time kaiju \
-v \
-z 5 \
-e 3 \
-m 11 \
-s 65 \
-t $KaDB/nodes.dmp \
-f $KaDB/kaiju_db_nr.fmi \
-i $KDOUT/${TEST}_bt2decon_R1.fastq.gz \
-j $KDOUT/${TEST}_bt2decon_R2.fastq.gz \
-o $KaOUT/${TEST}_kaiju



# sbatch $MAT/slurm_kaij_loop.sh $INES/Materials/samples $KDOUT $KaOUT $KaDB  10

# check for the first 6 samples done
for i in $( cat $INES/Materials/samples | head -9  | tail -3 );
do 
  KAI_C=$(grep -cE '^C' $KaOUT/${i}_kaiju )
  KAI_U=$(grep -cE '^U' $KaOUT/${i}_kaiju ) 
  KAI_TOT=$(echo "$KAI_C + $KAI_U" | bc)
  KAI_PC=$(echo "scale=2; ($KAI_C / $KAI_TOT)*100" | bc)
  echo "## kaiju sample processed: ${i} : total classified: ${KAI_PC}% (total: $KAI_TOT read pairs)  ------"
done


# kraken2 standard (can use bracken)
lk $INES/4__krak2

# kraken2 maxi (no bracken) - in jamie.dir
lk $MSDAT/r0936/4__krak2_maxi/*kraken_report

## all taxonomy together, *AT THE SPECIES LEVEL* - some issue solved by adding -l arg in kaiju2Table
kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -r species -c 10 -l domain,phylum,class,order,family,genus,species -o $KaOUT/fhi__redch__kaiju_speciesab.tsv $KaOUT/*_kaiju

less $KaOUT/total_kaiju__summary_species.tsv


# NUR  ================================================

# kraken2 using maxikraken
MXIKDB=$MSDAT/ref/maxikraken2_1903_140GB

module load kraken2/2.1.1

WRK=/data/Food/analysis/R0602_microsupport/jamie.fitzgerald/r0937
KDOUT=$WRK/3__knead
KROUT=$WRK/4__krak2
mkdir ${KROUT}_default
mkdir ${KROUT}_maxi
mkdir ${KROUT}_fonly

# test for recouping taxonomy
KR_threads=5

# test
time kraken2 --db $K2REFDIR \
${KDOUT}/${TEST}_bt2decon_R1.fastq.gz \
${KDOUT}/${TEST}_bt2decon_R2.fastq.gz \
--paired \
--confidence 0.5 \
--minimum-hit-groups 5 \
--minimum-base-quality 20 \
--report-zero-counts \
--gzip-compressed \
--threads $KR_threads \
--report ${KROUT}_default/${TEST}_test_kraken_report \
--unclassified-out ${KROUT}_default/${TEST}_test_kraken_unclass# \
--output ${KROUT}_default/${TEST}_test_kraken_output

# KrakenTools: convert to mpa stylee
~/bin/KrakenTools/kreport2mpa.py -r ${KROUT}_default/${TEST}_test_kraken_report -o ${KROUT}_default/${TEST}_kraken_mpa
# pipe mpa file to reformat
grep -h '|s_' ${KROUT}_default/${TEST}_kraken_mpa | cut -f 1 | sort | uniq | sed 's/|/\t/g' > ${KROUT}_default/krakenStnd_taxonomy.tsv

# spec k2 default as below, re-run with the 125 nt database in $MSDAT/ref/kraken2
K2REFDIR=$MSDAT/ref/kraken2
sbatch $MAT/slurm_krak2_loop.sh $NUR/Materials/samples $KDOUT ${KROUT}_default $K2REFDIR 10
sbatch $MAT/slurm_krak2_loop.sh $NUR/Materials/samples $KDOUT ${KROUT}_maxi $MXIKDB 10
sbatch $MAT/slurm_krak2-fonly_loop.sh $NUR/Materials/samples $NUR/2__filt ${KROUT}_fonly $MXIKDB 10

# bracken on kracken-def - dcheck you've right db version etc.
module load braken

BR_r=125
BR_l=S
BR_t=10   # counts! not threads

for i in $(cat $NUR/Materials/samples );
  do bracken -d $MSDAT/ref/kraken2 -i ${KROUT}_default/${i}_kraken_report -o ${KROUT}_default/${i}.bracken -r $BR_r -l $BR_l -t $BR_t ;
done                                                                                        

## combine reports to one monolith.tsv
combine_bracken_outputs.py --files ${KROUT}_default/*.bracken -o ${KROUT}_default/fhi_chickmi_krakenStnd_abundances.tsv


# kaiju standard   =======================================================================================

ls -lsh $MSDAT/r0937/4__kaiju

KaOUT=$MSDAT/r0937/4__kaiju
for i in $( cat $NUR/Materials/samples );
do 
  KAI_C=$(grep -cE '^C' $KaOUT/${i}_kaiju )
  KAI_U=$(grep -cE '^U' $KaOUT/${i}_kaiju ) 
  KAI_TOT=$(echo "$KAI_C + $KAI_U" | bc)
  KAI_PC=$(echo "scale=2; ($KAI_C / $KAI_TOT)*100" | bc)
  echo "## kaiju sample processed: ${i} : total classified: ${KAI_PC}% (total: $KAI_TOT read pairs)  ------"
done


# _________________________
# now convoluted ## ==========================
 ## ==========================


[...] GCA_003422985.1, identical, 
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/422/985/GCF_003422985.1_ASM342298v1, , , naGCF_003423005.1, PRJNA224116, SAMD00083667, BEDF00000000.1, na, 1280, 1280, Staphylococcus aureus, strain=M6K209, , latest, Scaffold, Major, Full, 2018/08/20, ASM342300v1, Medical Mycology ResearchSM2031333v1, , , nafailed: 1024 at /data/Food/analysis/R0602_microsupport/jamie.fitzgerald/ref/k_db_build_gh/download_bacteria.pl line 62, <IN> line 89762.
(bperl) [jamie.fitzgerald@compute05 maxikraken_nov22]$

 ## ==========================



TEST=_N5_Nuria_S296_L001

TEST=Q4T8


## new DB  ---------

## downaloaded...

## make a new bperl conda
conda create -m bperl
conda activate bperl
conda install -c bioconda perl-bioperl

# for git scripts
KT=$MSDAT/ref/k_db_build_gh

# for db
DB_NAME=$MSDAT/ref/maxikraken_nov22
mkdir $DB_NAME
cd $DB_NAME

#for genomes
mkdir $DB_NAME/k_genomes ; cd $DB_NAME/k_genomes

## dl from genomes folder
#archaea
perl $KT/download_archaea.pl &
#bacteria
perl $KT/download_bacteria.pl >> k2_bact.dl_overflow.stout &
#fungi
perl $KT/download_fungi.pl &
#protozoa
perl $KT/download_protozoa.pl &
#viral
perl $KT/download_viral.pl &
#human
perl $KT/download_human.pl &

# copy bt2 ref genomes to DB site (delete human duplicate) - note will need to be FNA - check decoy genome!
parallel cp {} $DB_NAME/ ::: /data/Food/primary/R0936_redcheese/ref/*fna.gz ; rm $MSDAT/ref/maxikraken_nov22/*GRCh38.p14*fna.gz
# unzip for matching
gzip -d $MSDAT/ref/maxikraken_nov22/*fna.gz

#check numbers
parallel "ls -lsh {} | wc -l" ::: $DB_NAME/k_genomes/*
# > [...] other single genomes - human, cow, goat, sheep, chicken, decoy
# > 97		# protozoa
# > 484		# fungi
# > 1339	# archaea	
# > 11740	# viruses
# > 29829	# bacto

# note that dl took 3 days, and didnt complete. get bact genome list:
# although note that that list is coming from refseq, while maxikrkaen is supposed to be intersect( refseq, genbank)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt


# nm all that, go parallel
#test 
find k_genomes/ -name '*.fna' -print | parallel -j 32 echo {} >> $KT/genomes.txt
less $KT/genomes.txt

module load kraken2/2.1.1
find k_genomes/ -name '*.fna' -print | parallel -j 32 bracken-build --add-to-library {} --db $DB_NAME
