
## try again! using shortbred to type the CARD database
	# watch out though, shortbred is nearly 10...
    
    # base
    CDB=~/ref/card
    MSDAT=$MSDAT    # sourced at login
    RAW=$MSDAT/raw/ms__paedcf_raw
    WRK=$MSDAT/ms__paedcf           #   no longer using ORPH flow, using main seqs (assembled for Kaiju)
    #
    # scripts etc
    PAED=~/ms__paedcf
    MAT=$PAED/Materials
    #
    # data processing etc
		# FILT=$MSDAT/ms__paedcf/2__filt     # picking it up here...
		# KDOUT=$WRK/3__knead
	  KaOUT=$WRK/4__kaiju
    # BTOUT=$WRK/6__bt
    
    # test
    TEST=SC145C1-NS_S86
    
  	mkdir $BTOUT   # assume all others precede this
	
    module load  parallel shortbred 
    module load  diamond  # sep out for conda complications
    module load  blast
    module load  usearch/6.1
    module load  muscle
    module load  cdhit
    module load  rapsearch
    
    
    ## get files
      mkdir ~/ref/card ; mkdir $BTOUT ; cd ~/ref/card 
      mkdir card_dat card_ont card_short
      # cd $DB/card_ont ; wget https://card.mcmaster.ca/download/5/ontology-v3.1.3.tar.bz2 ; tar -xvf ontology-v3.1.3.tar.bz2 &
      # cd $DB/card_dat ; wget https://card.mcmaster.ca/download/0/broadstreet-v3.1.3.tar.bz2 ; tar -xvf broadstreet-v3.1.3.tar.bz2 &
      cd $CDB/card_short ; wget https://github.com/biobakery/shortbred/releases/download/0.9.4/ShortBRED_CARD_2017_markers.faa.gz &
      cd $CDB
	  
	  less $CDB/card_short/ShortBRED_CARD_2017_markers.faa.gz
	  wc -c $CDB/card_short/ShortBRED_CARD_2017_markers.faa.gz   # 98 MB
	  wc -l $CDB/card_short/ShortBRED_CARD_2017_markers.faa.gz   # 200 markers! D:


    ## test files to index
      # https://github.com/biobakery/shortbred/
	  # http://huttenhower.sph.harvard.edu/shortbred
	  # 
	  
    ## test kit -----------------------------

		# shortbred-identify  - we skip this bit as CARD markers already provided for by the SHORTBRED bakers

		lk $KaOUT/${TEST}/${TEST}_R1_kneaddata_paired_1_stksort.fastq.gz
shortbred_quantify.py \
--markers $CDB/card_short/ShortBRED_CARD_2017_markers.faa.gz \
--wgs $KaOUT/$TEST/${TEST}__R1_kneaddata_paired_1_stksort.fastq.gz \
--results results.txt \
--marker_results $BTOUT/$TEST/${TEST}_marker.out \
--threads 5 \
--diamond /install/software/anaconda3.6.a/bin/diamond \
--tmp  $BTOUT/$TEST/${TEST}_tmp
	

    
    ## slurm ver    -------------------------
    echo '#!/bin/bash

				#SBATCH –-job-name=cardy_beats
				#SBATCH –-output=12xrphan_bt_card.txt

				#SBATCH –-ntasks=5
				#SBATCH –-time=3

				SAMPLE=$1
				IN=$2
				OUT=$3
				REFDIR=$4

				##   C A R D    b o w t i e 
				time bowtie2 --no-unal -p 5 --very-sensitive -x $DB/basic_bt2 -U $IN/${SAMPLE}/${SAMPLE}_R1R2Xrphan.fastq.gz -S $OUT/${SAMPLE}_card_bt2.sam

				
		' > $MAT/slurm_card_bt2.sh
    
    
    ## run on cohort
    ls $KDOUT/*/*_R1R2Xrphan.fastq.gz  | sed -e 's/.*\/\(.*\)_R1R2X..*/\1/g' | sort | uniq | parallel -j 8 sbatch $MAT/slurms/slurm_card_bt2.sh {} $KDOUT $BTOUT $DB
    
    # fast, because bad
    
    
    ## combine to feature table
    
        ... have no idea how to do this. 
        
      
        