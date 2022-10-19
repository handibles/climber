

##  B O W T I E   C A R D   D B   =========================================================

	## not that this is without merit, but
		# - a better choice is runnning this via shortBred, which will::
		# - - estimate abundances for you
		# - - search the protein sequences
		# - - perform a faster reduced sequence search

		
	## try BT2 some other time for some other task  ...


  ## recap vars
    
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
    FILT=$MSDAT/ms__paedcf/2__filt     # picking it up here...
    KDOUT=$WRK/3__knead
	KaOUT=$WRK/4__kaiju
    BTOUT=$WRK/6__bt
    
    # test
    TEST=SC145C1-NS_S86
    
	mkdir $BTOUT   # assume all others precede this
	
    module load parallel bowtie2 python3.7
    
    
    ## get files
      # mkdir ~/ref/card ; mkdir $BTOUT ; cd ~/ref/card 
      # mkdir card_dat card_ont
      # cd $CDB/card_ont ; wget https://card.mcmaster.ca/download/5/ontology-v3.1.3.tar.bz2 ; tar -xvf ontology-v3.1.3.tar.bz2 &
      # cd $CDB/card_dat ; wget https://card.mcmaster.ca/download/0/broadstreet-v3.1.3.tar.bz2 ; tar -xvf broadstreet-v3.1.3.tar.bz2 &
      # cd $CDB
      # 
    ## build index -  4am pain, but v fast
      # cd $CDB/card_dat
      # bowtie2-build nucleotide_fasta_protein_homolog_model.fasta,nucleotide_fasta_protein_knockout_model.fasta,nucleotide_fasta_protein_overexpression_model.fasta,nucleotide_fasta_protein_variant_model.fasta,nucleotide_fasta_rRNA_gene_variant_model.fasta,protein_fasta_protein_homolog_model.fasta,protein_fasta_protein_variant_model.fasta,protein_fasta_protein_overexpression_model.fasta,protein_fasta_protein_knockout_model.fasta $CDB/basic_bt2

	# check old SAM outputs! 
		ls /data/Food/primary/R0602_microsupport/ms__paedcf_orph/6__bt
		less -S /data/Food/primary/R0602_microsupport/ms__paedcf_orph/6__bt/${TEST}_card_bt2.sam
    
    ## test files to index
      # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner
      # bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | --sra-acc <acc> | b <bam>} -S [<sam>]
      # - --no-unal  :  drop unaligned reads
	  
	  # expect that we would like high specificity global alignments for this to be sensible - also, refs are short!
    lk $KaOUT/${TEST}/${TEST}_R1_kneaddata_paired_1_stksort.fastq.gz
	
    bowtie2 --no-unal --no-hd -p 5 -N 1 --very-sensitive -x $CDB/basic_bt2 -1 $KaOUT/${TEST}/${TEST}_R1_kneaddata_paired_1_stksort.fastq.gz -2 $KaOUT/${TEST}/${TEST}_R1_kneaddata_paired_2_stksort.fastq.gz -S $BTOUT/${KSAMPLE}_card_bt2_paired.sam
    
    bowtie2 --no-unal --no-hd -p 5 -N 1 --very-sensitive -x $CDB/basic_bt2 -U $KaOUT/${TEST}/${TEST}_R1_kneaddata_unmatched_12.fastq.gz -S $BTOUT/${KSAMPLE}_card_bt2_unmatched.sam
    
    
    ## slurm ver    
    echo '
#!/bin/bash
#SBATCH –-job-name=cardy_beats
#SBATCH –-output=12xrphan_bt_card.txt
#SBATCH –-ntasks=1
#SBATCH –-cpus-per-task=5
#SBATCH –-ntasks-per-node=1
#SBATCH –-time=3

SAMPLE=$1
IN=$2
OUT=$3
REFDIR=$4

mkdir $OUT/$SAMPLE

##   C A R D    b o w t i e 
bowtie2 --no-unal --no-hd -p 5 -N 1 --very-sensitive -x $CDB/basic_bt2 -1 $IN/${SAMPLE}/${SAMPLE}_R1_kneaddata_paired_1_stksort.fastq.gz -2 $IN/${SAMPLE}/${SAMPLE}_R1_kneaddata_paired_2_stksort.fastq.gz -S $OUT/${SAMPLE}_card_bt2_paired.sam

bowtie2 --no-unal --no-hd -p 5 -N 1 --very-sensitive -x $CDB/basic_bt2 -U $IN/${SAMPLE}/${SAMPLE}_R1_kneaddata_unmatched_12.fastq.gz -S $OUT/${SAMPLE}_card_bt2_unmatched.sam

cat $OUT/${SAMPLE}_card_bt2_*.sam > $OUT/${SAMPLE}_card_combo.tsv

' > $MAT/slurm_card_bt2.sh
    
    
    ## deal CARDS
    ls -d $KaOUT/*/  | sed -E 's/.*\/(.*)\/$/\1/g' | sort | uniq | parallel -j 22 sbatch $MAT/slurms/slurm_card_bt2.sh {} $KaOUT $BTOUT $CDB
    


## APPROACH  ====================================================

  ## ridiculously, STILL, have to then copy if to a local unix box, then USB it to here... 
    scp -i ~/.ssh/id_rsa_ucc $BTOUT/*tsv jamie@dunsh.ucc.ie:/home/jamie/
    ssh -i ~/.ssh/id_rsa_ucc jamie@dunsh.ucc.ie "scp -P 55555 ~/*tsv jfg@143.239.154.14:/home/jfg/Dropbox/teag_sync" & 
    

        
      
        