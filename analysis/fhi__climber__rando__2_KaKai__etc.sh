


## needs:
# graphical representations
# gneomes for database

# get kraken2 abundances merged.


## --------------

K2REFDIR=$MSDAT/ref   # 130mers
MXIKDB=$MSDAT/ref/maxikraken2_1903_140GB
KaDB=/data/databases/kaiju/nr
WRK=$MSDAT/r0936  # ines
#WRK=MSDAT/r0937  # nuria
KR_threads=5
KDOUT=$WRK/3__knead
KaOUT=$WRK/4__kaiju
KROUT=$MSDAT/r0936/4__krak2

TEST=Q4T8


## ================================================

module load kraken2/2.1.1
module load braken

# test/taxon
time kraken2 --db $K2REFDIR \
$KDOUT/${TEST}_bt2decon_R1.fastq.gz \
$KDOUT/${TEST}_bt2decon_R2.fastq.gz \
--paired \
--confidence 0.15 \
--report-zero-counts \
--minimum-hit-groups 3 \
--minimum-base-quality 10 \
--gzip-compressed \
--threads $KR_threads \
--report ${KROUT}_def15.3.10/${TEST}_kraken_report \
--output ${KROUT}_def15.3.10/${TEST}_kraken_output

# slurm
KR_threads=5
sbatch $MAT/slurm_krak2_loop.sh $INES/Materials/samples $MSDAT/r0936/3__knead $MSDAT/r0936/4__krak2_default $K2REFDIR $KR_threads

BR_r=130
BR_l=S
BR_t=10   # counts! not threads

for i in $(cat $INES/Materials/samples );
  do bracken -d $MSDAT/ref/ -i ${KROUT}_def15.3.10/${i}_kraken_report -o ${KROUT}_def15.3.10/${i}.bracken -r $BR_r -l $BR_l -t $BR_t ;
done                                                                                        

# bracken outputs
combine_bracken_outputs.py --files ${KROUT}_def15.3.10/*.bracken -o ${KROUT}_def15.3.10/fhi_chickmi_krakenStnd_abundances.def15.3.10.tsv

# taxonomy
~/bin/KrakenTools/kreport2mpa.py -r ${KROUT}_def15.3.10/${TEST}_test_kraken_report -o ${KROUT}_def15.3.10/${TEST}_kraken.def15.3.10_mpa 
grep -h '|s_' ${KROUT}_default/${TEST}_kraken.def15.3.10_mpa | cut -f 1 | sort | uniq  sed 's/|/\t/g' > ${KROUT}_default/fhi_redch_krakenStnd_taxonomy.def15.3.10.tsv


## ================================================

## redo kaiju for $TEST also
KaDB=/data/databases/kaiju/nr
ls -lsh $MSDAT/r0936/4__kaiju
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

kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -r species -c 10 -l domain,phylum,class,order,family,genus,species -o $KaOUT/fhi__redch__kaiju_speciesab.tsv $KaOUT/*_kaiju
less $KaOUT/total_kaiju__summary_species.tsv


## ================================================

TESTDIR=$SHAR/r0936/4__krak2_def15
for i in $(cat $INES/Materials/samples );
do
	echo $i $(grep -E 'unclassified$' $TESTDIR/${i}_kraken_report) ;
done > ines_krak_unclass_counter.txt

# check for the first 6 samples done
for i in $( cat $INES/Materials/samples );
do 
  KAI_C=$(grep -cE '^C' $KaOUT/${i}_kaiju )
  KAI_U=$(grep -cE '^U' $KaOUT/${i}_kaiju ) 
  KAI_TOT=$(echo "$KAI_C + $KAI_U" | bc)
  KAI_PC=$(echo "scale=2; ($KAI_C / $KAI_TOT)*100" | bc)
  echo "## kaiju sample processed: ${i} : total classified: ${KAI_PC}% (total: $KAI_TOT read pairs)  ------"
done > ines_kaiju_output.txt


## ================================================


## check proceeds
grep -E unclassified$ ${KROUT}_fonly/_N*report
grep -E unclassified$ ${KROUT}_default/_N*report
grep -E unclassified$ ${KROUT}_def15.3.10/_N*report
grep -E unclassified$ ${KROUT}_maxi/_N*report
grep -E unclassified$ ${KROUT}_maxi.15.3.10/_N*report



