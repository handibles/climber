    
    # vars
    DB=/data/databases
  	HGR38_BT=/data/databases/hostremoval/Homo_sapiens/Bowtie2/
    MSDAT=$MSDAT
    WRK=$MSDAT/ms__paedcf
    RAW=$MSDAT/raw/ms__paedcf_raw
    #
    # scripts etc
    PAED=~/ms__paedcf
    MAT=$PAED/Materials
    SLURM=$PAED/slurms
    #
    # data processing etc
    QC=$WRK/1__qc
    FQC=$QC/fastqc
    FILT=$WRK/2__filt
    KDOUT=$WRK/3__knead
    # KROUT=$WRK/4__krak2
    
    KaOUT=$WRK/4__kaiju
    # KaDB=/data/databases/kaiju_v1.7.4/nr_euk/nr_euk     # lets try nr-euk! 
    KaDB=/data/databases/kaiju/nr     # or not, its very big... 
    TEST=SC100C6-TS_S51                                 # tester - 3m20s


		module load parallel fastqc multiqc kaiju kneaddata



##   E M P T I E S   ===============================================

    # theres always something
    
    # losing samples at two steps - time-killed at kneadata ( >90min) 
    #                             - segfault for kaiju
    
    # need to grab unrepresented samples, and then pack back into pipes....
    

## grab missing   ===================================================================================================

    COMP=$MAT/compare ; mkdir $COMP

  ## raw materials
    ls $RAW/*gz |  sed -E 's/.*\/(.*)_L00.*/\1/g' | sort | uniq > $COMP/raw_temp
    wc -l $COMP/raw_temp

    
  ## files are filtered in trimmo :: $FILT
    ls $FILT/*r1R1* # 89
    ls $FILT/*r2R2* # 89
    ls $FILT/*r2R2* | sed -E 's/.*\/(.*)_r2R2.*/\1/g' | sort | uniq > $COMP/filt_r2R2_temp

  
  ## paired and unpaired reads reassembled into r1R1/r2R2 fastq.gz, then given to KD2
    ls $KDOUT/*/*paired_1.fastq.gz | wc -l      # ~~82 now all 88! 
    ls $KDOUT/*/*paired_2.fastq.gz | wc -l      # ~~81
    ls $KDOUT/*/*unmatched_1.fastq.gz | wc -l   # ~~81
    ls $KDOUT/*/*unmatched_2.fastq.gz | wc -l   # ~~81
    ls $KDOUT/*/*paired_2.fastq.gz | sed -E 's/.*\/(.*)_R{1,2}.*/\1/g' | sort | uniq > $COMP/kd2_paired_temp


  ## and then there are the pieces missing via kaiju / kraken
    ls $KaOUT/*/*tsv | sed -E 's/.*\/(.*)_kaiju.*/\1/g' | sort | uniq > $COMP/kai_summ_temp
    wc -l $COMP/kai_summ_temp    # 89

  ## or could simply have failed, but can't assume that the "empty" isn't from a previous run. 
  ## instead check the 0 size outputs
	lk $KaOUT/*/*.out | grep ' 0 Oct' | sed -E 's/.*\/(.*)_kaiju.*/\1/g' | sort | uniq > $COMP/kai_empty_temp
  
	
	
	
  ## oh wow not fun
  ## more reliable to simply look through all the outputs and check for empties (probably a more logical way)
			# # ls $KaOUT/*/*out

			# isfinished() {
				# FILE=$1
				# MINSI=1000
				# ACTUSI=$(wc -c <"$FILE")
				# if [ $ACTUSI -ge $MINSI ]; then
					# return $1
				# else
					# exit
				# fi
			# }
			# export -f isfinished





  ## find those truants
    FULLSAMP=$COMP/raw_temp
    # COMPLETED=$COMP/filt_r2R2_temp   # filt completed
    COMPLETED=$COMP/kd2_paired_temp     # KD2 missing 
    comm -23 $FULLSAMP $COMPLETED  > $COMP/comp_reload
    cat $COMP/comp_reload

    cat $COMP/comp_reload | parallel ls -lsh $RAW/{}*    


## and then...                                                  =======================================================
##              ... well ...
##                                    ...reassemble I guess ?

	## failed filter / QC

		# give the SLURM more tiiiime - from 90min to 180 ? 
		sed 's/--time=90:00/--time=180:00/g' $MAT/slurms/slurm_kd2.sh > $MAT/slurms/slurm_kd2_extended.sh


		## consider clearing the decks for these before setting up 
		# cat $COMP/comp_reload | rm -rf $KDOUT/{}
		## some jobs (1? 2?) still killed by time. refire. 
		
		# knead more  
			cat $COMP/comp_reload | parallel -j 9 sbatch $MAT/slurms/slurm_kd2_extended.sh {} $FILT $KDOUT $HGR38_BT $FQC/ms__paedcf_knead

		# christ, wait for the above to finish first ! 
			# 
			#     # cross those simian digits
				cat $COMP/comp_reload | parallel ls -lsh $KDOUT/{}/*
				cat $COMP/comp_reload | parallel -j 9 sbatch $MAT/slurms/slurm_kaiju.sh {} $KDOUT $KaOUT $KaDB


	## failed kaiju

				cat $COMP/kai_empty_temp | parallel -j 9 sbatch $MAT/slurms/slurm_kaiju_single.sh {} $KDOUT $KaOUT $KaDB


	#  --- simpler just to look them up and re-run
	
SAMPLE=LIBNEG_S89
SAMPLE=SC128C3-NS_S55   # still being processed
# SAMPLE=LIBNEG_S89
# SAMPLE=LIBNEG_S89
IN=$KDOUT
REFDIR=$KaDB
OUT=$KaOUT


time kaiju \
   -v \
   -z 5 \
   -e 3 \
   -m 11 \
   -s 75 \
   -t $KaDB/nodes.dmp \
   -f $KaDB/kaiju_db_nr.fmi \
   -i $KaOUT/${SAMPLE}/${SAMPLE}_R1_kneaddata_paired_1_stksort.fastq.gz \
   -j $KaOUT/${SAMPLE}/${SAMPLE}_R1_kneaddata_paired_2_stksort.fastq.gz \
   -o $KaOUT/${SAMPLE}/${SAMPLE}_kaiju__paired.out	>> $KaOUT/${SAMPLE}.stout

time kaiju \
   -v \
   -z 5 \
   -e 3 \
   -m 11 \
   -s 75 \
   -t $KaDB/nodes.dmp \
   -f $KaDB/kaiju_db_nr.fmi \
   -i $KaOUT/${SAMPLE}/${SAMPLE}_R1_kneaddata_unmatched_12.fastq.gz \
   -o $KaOUT/${SAMPLE}/${SAMPLE}_kaiju__unmatched.out
		
		

          
  
