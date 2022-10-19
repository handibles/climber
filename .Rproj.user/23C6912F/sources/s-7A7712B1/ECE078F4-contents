
## Ines' Red Cheese data

# check MD5

# check neg ctrl w.r.t. BAL samples  -  consider asking Raul



# ( ... get daata ...)

	# paralell --will cite ...
	
	# count files - 752?

	
## parameters for you

	# screen first!
	
	DB=/data/data/databases
	HGR38_BT=/data/databases/hostremoval/Homo_sapiens/Bowtie2/
	K2REFDIR=/data/databases/Kraken2_GTDB_r89_54k
	MSDAT=$MSDAT
	WRK=$MSDAT/fhi_redch
	RAW=$MSDAT/raw/fhi_redch_raw
	#
	# scripts etc
	REDCH=~/fhi_redch
	MAT=$REDCH/Materials
	SLURM=$REDCH/slurms
	#
	# data processing etc
	QC=$WRK/1__qc
	FQC=$QC/fastqc
	FILT=$WRK/2__filt
	KDOUT=$WRK/3__knead
	KROUT=$WRK/4__krak2

	module load parallel fastqc multiqc kraken2
	

	# data dir
	mkdir $MSDAT/fhi_redch $WRK/2__filt $WRK/3__knead $WRK/1__qc $QC/fastqc $FQC/fhi_redch_raw $FQC/fhi_redch_trimm $WRK/4__krak2
	# analysis dir
	mkdir ~/fhi_redch $REDCH/Materials $REDCH/slurms
	mkdir $FILT/fhi_redch_raw
	mkdir $FQC/fhi_redch_trimm $FQC/multi_trimm


## SLURM template   ---------------------------
	
		# slurm params https://bugs.schedmd.com/show_bug.cgi?id=2724
		# need to manually insert the # character
		# mem alloc doesn't work
		# https://slurm.schedmd.com/pdfs/summary.pdf
		
		# ## missing outputs
		# for i in $(ls $RAW/*gz | sed -e 's/.*\/\(.*\)_R._001.*/\1/g' | sort | uniq ) ; do ls $FILT/*$i*out ; done | less
		# # simply delete, and redo with more time :(
		

## FQC / MQC	-------------------------------

		# threads for FQC will be allocated 1CPU per file,
		# i.e. -t 10 on 3 FASTQ uses only 3 CPU
		
		module load fastqc
		module load multiqc
		# ls $RAW/SC145C1-NS_S86_L001*
		# time fastqc -t 4 $RAW/SC145C1-NS_S86_L001* -o $FQC/fhi_redch_raw
		# ls $RAW/*fastq.gz | sed -e 's/.*\/\(.*\)_L00.*/\1/g' | sort | uniq | wc -l
		
		
		echo '#!/bin/bash

				SBATCH –job-name=knead_fastq
				SBATCH –output=knead_fastq.txt

				SBATCH –ntasks=8
				SBATCH –time=2:35

				FQFILE=$1
				IN=$2
				OUT=$3
		
				##   F Q C  /   M Q C
				time fastqc -t 8 $IN/$FQFILE* -o $OUT

		' > $MAT/slurms/slurm_fqc.sh
		cat $MAT/slurms/slurm_fqc.sh
		
		sbatch $MAT/slurms/slurm_fqc.sh SC145C1-NS_S86_L00 $RAW $FQC/fhi_redch_raw
		ls $RAW/*fastq.gz | sed -e 's/.*\/\(.*\)_L00.*/\1/g' | sort | uniq \
		| parallel --will-cite sbatch $MAT/slurms/slurm_fqc.sh {}* $RAW $FQC/fhi_redch_raw

	# return
		mkdir $FQC/multi_raw
		time multiqc $FQC/fhi_redch_raw -o $FQC/multi_raw
		scp -r $FQC/multi_raw jamie@143.239.204.169:/home/jamie/teag_sync ; ssh jamie@143.239.204.169 "to_asp teag_sync"
		# then need to grab from yoda, or send yoda - unfinished fn to do so


##   a d a p t e r   s e q u e n c e s   d e t e c t e d  -   up to 75%

		## teagasc has a verrrrry well stocked reference file for FASTQ - using those refs for trimmo, also
		## via http://www.science.smith.edu/cmbs/wp-content/uploads/sites/36/2020/01/illumina-adapter-sequences-1000000002694-11.pdf
		## several of these match ref in trimmomatic package
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

		## test case
		# time fastqc -t 4 $RAW/SC145C1-NS_S86_* -o $FQC/fhi_redch_raw

		cat 
		# ls $MAT/fqc_trimm_ill_ref.fa 
		ls -lsh $RAW/SC145C1-NS_S86_L001_R*_001.fastq.gz
		fastqc $RAW/SC145C1-NS_S86_L001_R*_001.fastq.gz -o $FQC/
		# check just the illumina clip stage
		fastqc -t=2	$FILT -o $FQC/
		# from hpc
		scp $FQC/*html jamie@143.239.204.169:/home/jamie/teag_sync
		## send :: ssh to yo, pass to asp
		yo
		ssh jamie@143.239.204.169 "to_asp teag_sync"
 

##   t r i m m o m a t i c    =================
 
		## using reference fqc_trimmo_ill_ref.fa from above
 		# interesting case for Trimmo over cutadapt, as SE reads *need* adapter sequence, PE reads *do not* : https://dnatech.genomecenter.ucdavis.edu/faqs/when-should-i-trim-my-illumina-reads-and-how-should-i-do-it/
		# poss incorrect, but no less suggests trimmomatic for this case
		
		echo '#!/bin/bash

		#SBATCH –job-name=trimmo_raw
		#SBATCH –output=trimmo_raw.txt

		#SBATCH –ntasks=4
		#SBATCH –time=06:00
		
		SAMPLE=$1
		IN=$2
		OUT=$3
		REFDIR
		
		# trimmo, removed keppBothReads opt ; decent 90% retention
		trimmomatic PE $IN/${SAMPLE}_R1_001.fastq.gz $IN/${SAMPLE}_R2_001.fastq.gz $OUT/${SAMPLE}_R1_trimm.fastq.gz $OUT/${SAMPLE}_R1_trimm_unpaired.fastq.gz $OUT/${SAMPLE}_R2_trimm.fastq.gz $OUT/${SAMPLE}_R2_trimm_unpaired.fastq.gz HEADCROP:25 CROP:125 ILLUMINACLIP:$REFDIR/fqc_trimmo_ill_ref.fa:2:30:10:5 SLIDINGWINDOW:6:15 MINLEN:80 -threads 6 > $OUT/trimmo_${SAMPLE}.out #2>&1
		
		' > $MAT/slurms/slurm_trimmo.sh
		# assumes absolute path
		ls $RAW/L*fastq.gz | sed -e 's/.*\/\(.*\)_R._001.*/\1/g' | sort | uniq | parallel --will-cite -j 6 sbatch $MAT/slurms/slurm_trimmo.sh {} $RAW $FILT $MAT

	
		# check
		cat $FILT/*out | less
		ls $FILT/*trimm.fastq.gz | sed -e 's/.*\/\(.*\)_R._trimm.*/\1/g' | sort | uniq | parallel --will-cite -j 10 sbatch $MAT/slurms/slurm_fqc.sh {}* $FILT $FQC/fhi_redch_trimm
		rm $FQC/fhi_redch_trimm/*trimm_unpaired*      # bad reads from here eff up the output ; tho must be a better way
		time multiqc $FQC/fhi_redch_trimm -o $FQC/multi_trimm
		scp -r $FQC/multi_trimm jamie@143.239.204.169:/home/jamie/teag_sync ; ssh jamie@143.239.204.169 "to_asp teag_sync"

	
##   r e - c a t   s a m p l e s   ===========================
		
		# good luck trying to cram 8 files into a kneaddata pipe.
		# read numbers per file DONT match, but they also don't match in raws  .... which is odd?  ( its just because you sourced the parts via glom)
		# cat gz directly using cat : https://unix.stackexchange.com/questions/158941/how-to-combine-gunzipped-fastq-files#296551
		# cat $FILT/${KSAMPLE}*_R1_trimm.fastq.gz > $FILT/${KSAMPLE}_R1.fastq.gz
		ls $FILT/*trimm.fastq.gz | sed -e 's/.*\/\(.*\)_L00._.*/\1/g' | sort | uniq | parallel -j 10 "cat $FILT/{}*_R1_trimm.fastq.gz > $FILT/{}_R1.fastq.gz" 
		ls $FILT/*trimm.fastq.gz | sed -e 's/.*\/\(.*\)_L00._.*/\1/g' | sort | uniq | parallel -j 10 "cat $FILT/{}*_R2_trimm.fastq.gz > $FILT/{}_R2.fastq.gz"
		
	
##   k n e a d d a t a 2   =======================
	
		## kneaddata - using /data/databases/hostremoval/Homo_sapiens/
		# re-run trimmomatic disabled using --bypass-trim

		## trial define
		# HGR38_BT -  for bt2 
		KIN=$FILT
		KOUT=$KDOUT
		KSAMPLE=SC145C1-NS_S86 #_L00*
		REFDIR=$HGR38_BT
		QCDIR=$FQC/fhi_redch_knead	

		module load kneaddata trimmomatic   # needs trimmo, even if you don't
		mkdir $FQC/fhi_redch_knead
		
		# trial
		time kneaddata -t 4 -p 4 --bypass-trim -i $KIN/${KSAMPLE}_R1.fastq.gz -i $KIN/${KSAMPLE}_R2.fastq.gz -o $KOUT/${KSAMPLE}_2  -db $HGR38_BT --max-memory 45g  --remove-intermediate-output > $KOUT/knead_${KSAMPLE}_2.stout --trimmomatic /install/software/anaconda2.7.a/share/trimmomatic # ; gzip $KOUT/$KSAMPLE/*fastq # 2>&1
	
		echo '#!/bin/bash

				#SBATCH –-job-name=knead_fastq
				#SBATCH –-output=knead_fastq.txt

				#SBATCH –-ntasks=4
				#SBATCH –-time=40:00
				SAMPLE=$1
		
				IN=$2
				OUT=$3
				REFDIR=$4
				QCDIR=$5

				##   K N E A D D A T A
				# if youre going to change the output names, change for STOUT< GZIP, and FASTQC also
				time kneaddata -t 4 -p 4 --bypass-trim -i $IN/${SAMPLE}_R1.fastq.gz -i $IN/${SAMPLE}_R2.fastq.gz -o $OUT/${SAMPLE} -db $REFDIR --max-memory 45g > $OUT/knead_${SAMPLE}.stout 2>&1
				gzip $OUT/$SAMPLE/*fastq 

				## mkdir
				## mkbreakfast
				## within dir, do fastqc -o $QCDIR/ -t 4 $OUT/${SAMPLE}/*.fastq.gz
				
				
		' > $MAT/slurm_kd2.sh
		
		ls $FILT/*R[12].fastq.gz | sed -e 's/.*\/\(.*\)_R..*/\1/g' | sort | uniq | parallel -j 6 sbatch $MAT/slurms/slurm_kd2.sh {} $FILT $KDOUT $HGR38_BT $FQC/fhi_redch_knead

		ls $FILT/*R[12].fastq.gz | sed -e 's/.*\/\(.*\)_R..*/\1/g' | sort | uniq | parallel ls $FILT/{}_R1.fastq.gz
	

##   d e a l   w i t h   h a l f  -  r u n

		ls $FQC/fhi_redch_knead | sed -e 's/\(.*\)_R..*/\1/g' | sort | uniq	
		
		mkdir $KDOUT/kd_stout ; mv $KDOUT/*out $KDOUT/kd_stout/
		
		mkdir $FQC/multi_knead
		multi -o $FQC/multi_knead $FQC/fhi_redch_knead
		scp -r $FQC/multi_knead jamie@143.239.204.169:/home/jamie/teag_sync ; ssh jamie@143.239.204.169 "to_asp teag_sync"

	# compare to find stragglers
		find $KOUT/*gz | sed -r 's/.*[\/_]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq > $DAN/Materials/temp2
		find $KRK/*report | sed -r 's/.*[\/_]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq > $DAN/Materials/temp3_k
		comm -23 $DAN/Materials/temp2 $DAN/Materials/temp3_k > $DAN/Materials/krk_offenders.txt  
		# include switch for outputting zero-reads - for tabulation




##   backup the important pieces... 2GB  --------------------------------------------------------------------------------

		ls $KDOUT/*/*_kneaddata_paired_[12].fastq.gz
		# parallel -j 10 scp {} jamie@143.239.204.169:/home/jamie/teag_sync ::: $KDOUT/*/*_kneaddata_paired_[12].fastq.gz
		parallel -j 10 scp {} jamie@143.239.204.169:/mnt/yabba/yoda_backup/home/ClaessonFTP/jamie/ ::: $KDOUT/*/*_kneaddata_paired_[12].fastq.gz



##   r e - r e - c a t   s a m p l e s 

		# kraken would like to run all the necessary data in one go - so append paired data with xxxxx / NNNNNNNN spacers, then cat with unpaired data, and run a single fastq.
		# ls $FILT/*R[12].fastq.gz | sed -e 's/.*\/\(.*\)_R..*/\1/g' | sort | uniq | parallel -j 20 
		usearch -fastq_join $KDOUT/$SAMPLE/${SAMPLE}_R1_kneaddata_paired_R1.fastq.gz -reverse $KDOUT/$SAMPLE/${SAMPLE}_R1_paired_R2.fastq.gz -fastaout $KDOUT/${SAMPLE}_kneaddata_R1R2.fastq.gz



##    k r a k e n    2   ===========================================================================================

		# test
		K2SAMPLE=SC145C1-NS_S86	
		K2IN=$KDOUT
		K2OUT=$KROUT
		K2REFDIR=/data/databases/Kraken2_GTDB_r89_54k

	# recall - first run requires a load to shared mem step, which takes time  -  to do this, possibly omit the memory mapping step. GOing to omit entirely as mem not an issue, yet
		cat $K2IN/${K2SAMPLE}/${K2SAMPLE}_R1_kneaddata_Homo*.fastq.gz > $K2IN/${K2SAMPLE}/${K2SAMPLE}_contam.fastq.gz &
		cat $K2IN/${K2SAMPLE}/${K2SAMPLE}_R1_kneaddata_unmatched_[12].fastq.gz > $K2IN/${K2SAMPLE}/${K2SAMPLE}_orphan.fastq.gz &
		# THREE MINUTES! only 18% classified though.  | 0.5 conf: 2m, 2% class | 0.15 cutoff - 3m, 13% ID'd
		time kraken2 --db $K2REFDIR  $K2IN/${K2SAMPLE}/${K2SAMPLE}_R1_kneaddata_paired_1.fastq.gz $K2IN/${K2SAMPLE}/${K2SAMPLE}_R1_kneaddata_paired_2.fastq.gz --paired --confidence 0.15 --gzip-compressed --threads 5 --report-zero-counts --report $K2OUT/${K2SAMPLE}_test_kraken_report# --unclassified-out $K2OUT/${K2SAMPLE}_test_kraken_unclass# --output $K2OUT/${K2SAMPLE}_test_kraken_output#
		# TWO-POINT-FIVE MINUTES! 27% identified  | conf 0.5: 1m30, 2% identifiedi | 3 MIN 30 - 15% classified
		time kraken2 --db $K2REFDIR  $K2IN/${K2SAMPLE}/${K2SAMPLE}_orphan.fastq.gz --gzip-compressed --confidence 0.15 --threads 5 --report-zero-counts --report $K2OUT/${K2SAMPLE}_test_orphan_kraken_report --unclassified-out $K2OUT/${K2SAMPLE}_test_orphan_kraken_unclass --output $K2OUT/${K2SAMPLE}_test_orphan_kraken_output
		# THREE MIN FORTY! 12% sequences identified, highly worrying.
		time kraken2 --db $K2REFDIR  $K2IN/${K2SAMPLE}/${K2SAMPLE}_contam.fastq.gz --gzip-compressed --confidence 0.15 --threads 5 --report-zero-counts --report $K2OUT/${K2SAMPLE}_test_contam_kraken_report --unclassified-out $K2OUT/${K2SAMPLE}_test_contam_kraken_unclass --output $K2OUT/${K2SAMPLE}_test_contam_kraken_output

		# need to specify id cutoff

		# for ROUND 1, paired-mode only is fine, as some disagreement about handling btw forum and manual - possibly has been updated in latter.
		# can remember why we included the zero counts...but we did, intentionally

		echo '#!/bin/bash

				#SBATCH –-job-name=knead_fastq
				#SBATCH –-output=knead_fastq.txt

				#SBATCH –-ntasks=5
				#SBATCH –-time=40:00

				SAMPLE=$1		
				IN=$2
				OUT=$3
				REFDIR=$4

				time kraken2 --db $REFDIR  $K2IN/${K2SAMPLE}/${K2SAMPLE}_R1_kneaddata_paired_1.fastq.gz $K2IN/${K2SAMPLE}/${K2SAMPLE}_R1_kneaddata_paired_2.fastq.gz --paired --confidence 0.15 --gzip-compressed --threads 5 --report-zero-counts --use-mpa-style --report $K2OUT/${K2SAMPLE}_kraken_report# --unclassified-out $K2OUT/${K2SAMPLE}_kraken_unclass# --output $K2OUT/${K2SAMPLE}_kraken_output#

				##
					## here, should cat and grab the orphan output 
				##
				
				scp $OUT/${SAMPLE}* jamie@143.239.204.169:/mnt/yabba/yoda_backup/home/ClaessonFTP/jamie/

				## within dir, do fastqc -o $QCDIR/ -t 4 $OUT/${SAMPLE}/*.fastq.gz
				
				
		' > $MAT/slurm_krak2.sh
		
		ls $FILT/*R[12].fastq.gz | sed -e 's/.*\/\(.*\)_R..*/\1/g' | sort | uniq | head -20 | tail -1 | parallel -j 8 sbatch $MAT/slurms/slurm_krak2.sh {} $KDOUT $KROUT $K2REFDIR
		ls $FILT/*R[12].fastq.gz | sed -e 's/.*\/\(.*\)_R..*/\1/g' | sort | uniq | head -20 | tail -1 | parallel -j 8 $MAT/slurms/bash_krak2.sh {} $KDOUT $KROUT $K2REFDIR

  
  	# lifeboat protocol
  		# K2SAMPLE=SC145C1-NS_S86	
  		K2IN=/claesson/jamie/fhi_redch
  		K2OUT=/claesson/jamie/fhi_redch/krak2
  		K2REFDIR=/claesson/jamie/ref/kraken2_stnd


      # << < < < need to add step here to start with one sample, removing --mem-map, and updating dirs/
			## remove --mpa-style as interferes with Bracken
			ls /claesson/jamie/fhi_redch/*[12].fastq.gz | sed -e 's/.*\/\(.*\)_R..*/\1/g' | sort | uniq | parallel -j 6 time kraken2 --db $K2REFDIR  $K2IN/{}_R1_kneaddata_paired_1.fastq.gz $K2IN/{}_R1_kneaddata_paired_2.fastq.gz --memory-mapping --paired --confidence 0.15 --gzip-compressed --threads 5 --report-zero-counts --report $K2OUT/{}_kraken_report# --unclassified-out $K2OUT/{}_kraken_unclass# --output $K2OUT/{}_kraken_output# >> ~/para_krak_fhi_redch.stout 2>&1
      
    # this'll put out all the seq too. if you only want the report, send your output/unclassified ..away 
    # try catch output, or anything else - because where does it go?
    
    ## pre-run very possibly not necessary with --memorymapping opt



##   B R A C K E N   B E A T I N G       ## ===================

      brthreads=10
      b=$DAN/3_b2  ; mkdir $b
      binst=~/BioApps
  
  		# K2SAMPLE=SC145C1-NS_S86	
  		WRK=/claesson/jamie/fhi_redch
  		K2IN=/claesson/jamie/fhi_redch/knead
  		K2OUT=/claesson/jamie/fhi_redch/krak2
  		K2REFDIR=/claesson/jamie/ref/kraken2_stnd

      BRK=$WRK/brack ; mkdir $BRK
      
    ## build data base first, using length 80bp 35kmers, takes an hour on 25 threads

    ## fast as immediate  -  run bracken on the kraken_report  -  note no threads arg!
          # -r = read lengths (taken as the starting read length, we're at about 125 by now...)
          # -l = Level, def=S ; (K, P, C, O, F, G, S )
          # -t = threshold, def = 10 ; level below which a taxo-branch on the kraken tree will be dropped from further assignment of abundances 
      ls $K2OUT/* | sed -r 's/.*\/(.*)_kraken.*/\1/g' | sort | uniq | \
        parallel -j 3  \
          python ~/BioApps/Bracken/src/est_abundance.py -i $K2OUT/*{}*report -k $K2REFDIR/database80mers.kmer_distrib -l S -t 10 -o $BRK/{}.bracken \
            > $BRK/bracken_output_82.stout
          # bracken -d $DB/kraken2_stnd -i $KRK/*{}_*report -o $b/{}_bracken_report.txt -r 150 -l S -t 10 > $b/bracken.stout

    ## combine reports to one monolith.tsv
      python2 ~/BioApps/Bracken/analysis_scripts/combine_bracken_outputs.py --files $BRK/*bracken -o fhi_redch_krakenStnd_bracken_combo
      
    ## convert to MPA-style report outputs, in order to get taxonomic tree (again, surely a better way than this)
      parallel python2 ~/BioApps/Bracken/src/kreport2mpa.py -r {} -o $BRK/{.}.mpa ::: $K2OUT/*report* > $BRK/fhi_redch_bragglers


  # # catch done/undone  
  #   find $KOUT/*gz | sed -r 's/.*[\/_]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq > $DAN/Materials/temp2
  #   find $KRK/*report | sed -r 's/.*[\/_]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq > $DAN/Materials/temp3_k
  #   comm -23 $DAN/Materials/temp2 $DAN/Materials/temp3_k > $DAN/Materials/krk_offenders.txt  
  #   # include switch for outputting zero-reads - for tabulation
  # 

##   T A X O N O M Y   -   glom from mpa output!   ## ========================================
      
    ## own funciton gettax_kraken : hacky, ++slow, and only gets species name. Avoid.  
    
    ## get all the S level assignments; trim counts ; sort; uniq
      grep -h '|s_' $BRK/*mpa | cut -f 1 > $BRK/bracken_82_taxa.temp
      sort $BRK/bracken_82_taxa.temp | uniq > $BRK/bracken_82_taxa.temp2
      sed 's/|/\t/g' $BRK/bracken_82_taxa.temp2 > $BRK/fhi_redch_krakenStnd_bracken_combo_taxa_82
  
    ## send home
      scp ./bracken_545_*[combo,taxa] $DSK

    # consider total taxonomy, but prob doesn't tell us anything re: missing levels
      cat $b/bracken_reports_545/*mpa | cut -f 1 > $b/bracken_545_total_tax.temp
      sort  $b/bracken_545_total_tax.temp | uniq > $b/bracken_545_total_tax.temp2
