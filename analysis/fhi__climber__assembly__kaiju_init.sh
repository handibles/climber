
## setup kaiju  ==========================================
    
conda create -n kaijamie
conda install -c bioconda kaiju
conda activate kaijamie

KaDB=/data/databases/kaiju/nr     # nr-euk very big..
Ka_THREADS=10
KaOUT=$WRK/4__kaiju
mkdir $KaOUT
# ulimit -s unlimited    # nasa mem trick


# TEST=Q4T8
TEST=_N5_Nuria_S296_L001


## test kaiju  ==========================================
    
# explicitly default
  #-e is allowed mismatches   - def = 3
  #-m is match length     - def = 11 
  #-s is match score       - def = 65

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

less -S $KaOUT/${TEST}_kaiju

KAI_C=$(grep -cE '^C' $KaOUT/${TEST}_kaiju )
KAI_U=$(grep -cE '^U' $KaOUT/${TEST}_kaiju )
echo "scale=3; $KAI_C / $KAI_U" | bc


## slurm-loop  kaiju   =======================================

echo '#!/bin/bash

#SBATCH --job-name=kaij_loop
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
sbatch $MAT/slurm_kaij_loop.sh  $NUR/Materials/samples $KDOUT $KaOUT $KaDB $Ka_THREADS ;


## count amounts
for i in $( cat $NUR/Materials/samples );
do 
  KAI_C=$(grep -cE '^C' $KaOUT/${i}_kaiju )
  KAI_U=$(grep -cE '^U' $KaOUT/${i}_kaiju ) 
  KAI_TOT=$(echo "$KAI_C + $KAI_U" | bc)
  KAI_PC=$(echo "scale=2; ($KAI_C / $KAI_TOT)*100" | bc)
  echo "## kaiju sample processed: ${i} : total classified: ${KAI_PC}% (total: $KAI_TOT read pairs)  ------"
done # >> delete_temp


## all taxonomy together - some issue solved by adding -l arg in kaiju2Table

kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -r genus -c 10 -l domain,phylum,class,order,family,genus,species -o $KaOUT/total_kaiju__summary_species.tsv $KaOUT/*_kaiju

# step off to R!



##   kaiju   ======================   ======================   ======================   ======================   ======================

    #   # https://kaiju.binf.ku.dk
    #   # https://github.com/bioinformatics-centre/kaiju     # more likley    
    #   # https://teaching.healthtech.dtu.dk/22126/index.php/Kaiju_exercise
    #   
    #   ## copy the DB to the existing /dev/shm/ (or mount & cp )
    #   ## alt : mount a memory device, and point DB to that : https://www.howtoforge.com/storing-files-directories-in-memory-with-tmpfs
    #   # ls $KaDB/*{nodes,names,fmi}*
    #   # TMP=/dev/shm
    #   # time cp $KaDB/*{nodes,names,fmi}* $TMP/ & 
    #   
    #   ## likely initial DB mapping step also, substantial if nr_euk
    #   ## who knows what defaults actually mean:
    #   
    #   ## tried stating memory needs (--mem) and changing to nr; tried stating ntask=1, ncpus = 5
    #   
    #     # nano $MAT/slurms/slurm_kaiju.sh
    #     # SAMPLE=$TEST
    #     # IN=$KDOUT
    #     # REFDIR=$TMP
    #     # OUT=$KaOUT
    # 
    #   ## removed sbatch arg for --mem, specific nodes
    # 
    # 
    # 
  # ## taxonomy   ==================================================
    # 
    #           #     ## should probably put this within SLURM - gets correct NA for blank ranks, but lost at summary step
    #           #     # parallel -j 22 kaiju-addTaxonNames -t $KaDB/nodes.dmp -n $KaDB/names.dmp -i {} -o {}.tax -r superkingdom,phylum,class,order,family,genus,species ::: $KaOUT/[L,S]*/*total 
    #           #   
    #           #   ## needs -p  ; still doesn't regularise the taxonomic output. 
    #           #   ## however, note also that this will drop a load of sequences if unclear where they belong  -  maybe don't go with species?
    #           #     parallel -j 22 kaiju2table -t $KaDB/nodes.dmp -n $KaDB/names.dmp -p -r species -o {}.tsv {} ::: $KaOUT/*/*tax
    #           #     cat $KaOUT/*/*tax.tsv | sed -E 's/.*\/(.*_kaiju__pup.total.*$)/\1/' | sed -E 's/_kaiju__pup.total.tax//g' > $KaOUT/ms__paedcf_kaiju.master.tax
    #           # 
    #           # 
    #           #     cat $KaOUT/*/*tax.tsv | sed -E 's/.*\/(.*_kaiju__pup.total).*/\1/' | head
    # 
    # 
  # ## APPROACH  ====================================================
    # 
    #   ## ridiculously, have to then copy if to a local unix box, then USB it to here... 
    #     scp -i ~/.ssh/id_rsa_ucc $KaOUT/*master* jamie@dunsh.ucc.ie:/home/jamie/
    #     ssh -i ~/.ssh/id_rsa_ucc jamie@dunsh.ucc.ie "scp -P 55555 ~/*master* jfg@143.239.154.14:/home/jfg/Dropbox/teag_sync"
    #     
    # 
    #     # then from dunshire to asp
  

