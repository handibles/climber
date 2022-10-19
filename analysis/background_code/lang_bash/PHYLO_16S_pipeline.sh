#!/bin/bash


#   <  !  >   T R A C K   Y O U R   R E A D S   <  !  >



#   v1.2.1  -  JAN_20RUN :: revisit to ensure primers removed 
#   v1.3.0  -  generalised github sb_4.11 version - project PHYLO, user USER
#              - keep in mind that dir structure etc still jamie/jfg and asp/biolinux specific
#              - keep in mind that IPs also visible
#              - consider renaming project STC



# ## =   ==  === =   ==  === =   ==  === =   ==  === =   ==  === =   ==  === =   ==  === =   ==  === =   ==

    ## 0) P R E P

    # need to change Rscript library path, or keep re-setting location using
          # # https://stackoverflow.com/questions/27673000/rscript-there-is-no-package-called
          # .libPaths(c("/home/jfg/R/x86_64-pc-linux-gnu-library/3.5", .libPaths()))

    # F I X:  set options to script for source/path/dir


# ## ============================================================================================================

  # scp jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiom__PHYLO/analyses/* ~/PHYLO/  ;  dos2unix ~/PHYLO/*


                # # # # # # # #
                # # # # # # #
                # # # # # # # #
                # # # HERE! #
                # # # # # # # #
                # # # # # # #
                # # # # # # # #



    if [ -d ~/Materials ] ; then  mv ~/Materials ~/completed_PHYLO_$(date "+%b_%d_%Y") ;  fi
    if [ -d /data/jamie/PHYLO/Materials ] ; then  mv /data/jamie/PHYLO/Materials /data/jamie/PHYLO/completed_PHYLO_$(date "+%b_%d_%Y") ; fi


## 0   -   ====================================================================

    # initial data commit - in this situation, no need to duplicate in your folder like a dad: overwrite with processed files

    if [ -d /data/jamie/PHYLO/Materials/0_raw_reads ]
    then
      echo 'r a w   r e a d s   a l r e a d y   i n   / d a t a / j a m i e   d i r .' \n  'C h e c k   f o r m a t   a n d   r u n #'
    else
      echo 'c p   r a w   r e a d s   t o   / d a t a / j a m i e   d i r'
      mkdir /data/jamie/PHYLO/Materials
      mkdir /data/jamie/PHYLO/Materials/0_raw_reads
 	    cp /data/jamie/jfg_rawdata/PHYLO_flotilla_1/* /data/jamie/PHYLO/Materials/0_raw_reads    ## run set by double wildcard: run*/*gz
    fi




  	# standardise SAMPLENAMES - sampl-e-n4m3.R1.fastsq  -  need to rationalise this for PHYLO and generally
      rename s/_L001__001_forward/.R1/g /data/jamie/PHYLO/Materials/0_raw_reads/*fastq -v
      rename s/_L001__001_reverse/.R2/g /data/jamie/PHYLO/Materials/0_raw_reads/*fastq -v
  	# ha! need better solution to sed renaming the path - underscore in dir names
    	cd /data/jamie/PHYLO/Materials/0_raw_reads
      rename s/_/-/g ./*fastq -v

    ## make a sample_run_list.txt of sample names and run associated - in this case run is always 1.
      # sed :  _line based_ editor , matching between lines problematic. Use tr to remove newlines once added, or awk to match across lines.
      # see also https://stackoverflow.com/questions/27510462/how-can-i-remove-double-line-breaks-with-sed
      mkdir ~/Materials
      cd ~/Materials

      echo -e 'SampleID\tRun' > sample_run_list.txt
      ls /data/jamie/PHYLO/Materials/0_raw_reads/ | sed s/.*R2\\.fastq//g >> sample_run_list.txt        # or (...| xargs -n1 basename |  cut -f 1 -d ".")
      sed s/\\.R1\\.fastq/\\t1/g sample_run_list.txt -i
      sed '/^$/d' sample_run_list.txt -i     #  ^ - line start ; $ - line end ; d flag = delete
      # scp ~/Materials/sample_run_list.txt jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiom__PHYLO/output


	#  =============================================================================


## 0.1  ::  Q U E R Y   L E N G T H S

 # file, total reads, read length*count  ::  see https://www.biostars.org/p/72433/ for awk trick

  cd /data/jamie/PHYLO/Materials/0_raw_reads/
  for i in $(ls /data/jamie/PHYLO/Materials/0_raw_reads/*.fastq | cut -f 7 -d "/");    # consider -Sr for smallest first   ;    < !!! >  beware -f # if path-level changes
  do echo $i ;                       # name
     cat $i | grep -c '^+$' ;        # read count ::   ^ , $ = BOL/EOL
     cat $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c ; # unique read lengths, per 4-line piece
  done > ~/Materials/0_read_profile.txt

  #scp ~/Materials/0_read_profile.txt jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiom__PHYLO/read_ouput



	# ==============================================================================================================================================




## 1) F A S T Q C ::    ==========================================================================================================================
   # https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/

  # # uncommment if actually using
  # cat /data/jamie/PHYLO/Materials/0_raw_reads/*R1.fastq > /data/jamie/PHYLO/Materials/0_raw_reads/PHYLO_raw_run1_R1.fastq
  # cat /data/jamie/PHYLO/Materials/0_raw_reads/*R2.fastq > /data/jamie/PHYLO/Materials/0_raw_reads/PHYLO_raw_run1_R2.fastq
  # fastqc --threads 18 /data/jamie/PHYLO/Materials/0_raw_reads/PHYLO_raw_run1_R1.fastq
  # fastqc --threads 18 /data/jamie/PHYLO/Materials/0_raw_reads/PHYLO_raw_run1_R2.fastq
  # mkdir /data/jamie/PHYLO/Materials/1_fastqc_out ; mv /data/jamie/PHYLO/Materials/0_raw_reads/*zip /data/jamie/PHYLO/Materials/0_raw_reads/*html /data/jamie/PHYLO/Materials/1_fastqc_out
  # rm /data/jamie/PHYLO/Materials/0_raw_reads/PHYLO_raw*
  # # scp -r /data/jamie/PHYLO/Materials/1_fastqc* jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiom__PHYLO/output



## 2) T R I M M O M A T I C ::    ================================================================================================================

  # http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

  # PE filtering = keeping DADA2 happy at cost of sensitivity

  mkdir /data/jamie/PHYLO/Materials/2_trimmed ; mkdir /data/jamie/PHYLO/Materials/2_trimmed/2.1_lost
  cd /data/jamie/PHYLO/Materials/0_raw_reads
  for i in $(ls *1.fastq | cut -f 1 -d ".") ;
    do java -jar /home/jamie/bin/trimmomatic.jar PE $i.R1.fastq $i.R2.fastq ../2_trimmed/trimmed_$i.R1.fastq ../2_trimmed/2.1_lost/lost_$i.R1.fastq ../2_trimmed/trimmed_$i.R2.fastq ../2_trimmed/2.1_lost/lost_$i.R2.fastq CROP:280 HEADCROP:19 LEADING:24 TRAILING:24 MINLEN:150 ;
  done

      # Check output - make a fn?
      cat /data/jamie/PHYLO/Materials/2_trimmed/*R1.fastq > /data/jamie/PHYLO/Materials/2_trimmed/PHYLO_trimm_R1.fastq
      cat /data/jamie/PHYLO/Materials/2_trimmed/*R2.fastq > /data/jamie/PHYLO/Materials/2_trimmed/PHYLO_trimm_R2.fastq
      fastqc --threads 26 /data/jamie/PHYLO/Materials/2_trimmed/PHYLO_trimm_R1.fastq
      fastqc --threads 26 /data/jamie/PHYLO/Materials/2_trimmed/PHYLO_trimm_R2.fastq
      mkdir /data/jamie/PHYLO/Materials/1_fastqc_out
      mv /data/jamie/PHYLO/Materials/2_trimmed/*fastqc* /data/jamie/PHYLO/Materials/1_fastqc_out ; rm /data/jamie/PHYLO/Materials/2_trimmed/PHYLO_trimm*
    #    scp -r /data/jamie/PHYLO/Materials/1_fastqc_out jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiom__PHYLO/output  #handy but needs key entry so dont script



## 3)  D A D A   2  ::  fastqPairedFilter =======================================================================================================

  ##  ensure that TRIMMOMATIC OUTPUT is reflected in truncLen parameters here - see FASTQCout: Sequence Length Distribution

  mkdir /data/jamie/PHYLO/Materials/3_d2
  cd /data/jamie/PHYLO/Materials/2_trimmed
  for i in $( ls *.R1.fastq  | cut -f 1 -d ".") ;
    do echo ~/R/bin/Rscript ~/PHYLO/analyses/3.8_DADA2.quality.R --inpath /data/jamie/PHYLO/Materials/2_trimmed --outpath /data/jamie/PHYLO/Materials/3_d2  --samplename $i  ;
  done > ~/Materials/step_3_fqPFilter.R

  parallel -j 16 < ~/Materials/step_3_fqPFilter.R   2>&1 | tee -a ~/Materials/3_d2_fqPFilter.log


      # check:
      cat /data/jamie/PHYLO/Materials/3_d2/*R1.fastq > /data/jamie/PHYLO/Materials/3_d2/PHYLO_dada_R1.fastq
      cat /data/jamie/PHYLO/Materials/3_d2/*R2.fastq > /data/jamie/PHYLO/Materials/3_d2/PHYLO_dada_R2.fastq
      fastqc --threads 18 /data/jamie/PHYLO/Materials/3_d2/PHYLO_dada*
      rm /data/jamie/PHYLO/Materials/3_d2/PHYLO_dada_R*.fastq
      mv /data/jamie/PHYLO/Materials/3_d2/PHYLO_dada* /data/jamie/PHYLO/Materials/1_fastqc_out
      # scp -r /data/jamie/PHYLO/Materials/1_fastqc_out jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiom__PHYLO/output/
      


## 4)  D A D A   2  ::   Error PER-RUN (4.4)    ==================================================================================================
  # need ~10^9 bp to accurately assess erroor profile for DADA -  IDEALLY take enough _random_ samples, in parallel, to reach 10^9 reads
  # see BNEJNEB approach (ITS?) 

  mkdir /data/jamie/PHYLO/Materials/4_error
  cd /data/jamie/PHYLO/Materials/3_d2
  
  ~/R/bin/Rscript ~/PHYLO/analyses/4.4_DADA2.error_pool.R --indir /data/jamie/PHYLO/Materials/3_d2 --outdir /data/jamie/PHYLO/Materials/4_error --samplelist ~/Materials/sample_run_list.txt --run 1



## 5)  D A D A   2  ::    Average Error :: Deprecated   ========================================================================================== 
  #
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.
  # DOUBLE CHECK that this is a good idea.



## 6)  D A D A   2  ::    DADA2 Final.         =================================================================================================== 
  # Dereplication, application of error profiles (step 4), chimeras checked, then pairs are merged to ASVs / sample.

  mkdir /data/jamie/PHYLO/Materials/5_final
  cd /data/jamie/PHYLO/Materials/3_d2
  
  for i in $(ls *.R1.fastq | cut -f 2 -d ".");
    do echo ~/R/bin/Rscript ~/PHYLO/analyses/6.2_DADA2.final.R --forwardReadpath /data/jamie/PHYLO/Materials/3_d2/d2.$i.R1.fastq --reverseReadpath /data/jamie/PHYLO/Materials/3_d2/d2.$i.R2.fastq --samplename $i --outdir /data/jamie/PHYLO/Materials/5_final --samplelist ~/Materials/sample_run_list.txt >> ~/Materials/step_6_final.R ;
  done
  parallel -j 14 <  ~/Materials/step_6_final.R  2>&1 | tee -a ~/Materials/6_d2_final.log



## 7)  D A D A   2  ::   R D S   A G G:
  # RDS files for each sample are length-filtered and aggregated into an overall sequence table and used to 'speculate' taxonomies 
  cd /data/jamie/PHYLO/Materials/5_final
  mkdir /data/jamie/PHYLO/Materials/6_PHYLO_Seqtab 
  mkdir /data/jamie/PHYLO/Materials/5_final/5_final_rdata ; mv /data/jamie/PHYLO/Materials/5_final/*RData /data/jamie/PHYLO/Materials/5_final/5_final_rdata
 
  ## needs to be run for each run, i.e. :: for i in $RUN_GLOBVAR  -   if script nukes all sample accessions, not handling runX_samplID text properly
  ~/R/bin/Rscript ~/PHYLO/analyses/7.3_DADA2.aggRDS.R --study 'PHYLO' --strand 'top' --run 1 --thresh 200 --datain '/data/jamie/PHYLO/Materials/5_final/' --taxref '/data/jamie/ref/SILVA_SSU_r132_March2018.RDS' --dataout '/data/jamie/PHYLO/Materials/6_PHYLO_Seqtab' 2>&1 | tee -a ~/Materials/7_d2_aggregate.log

  # scp -r /data/jamie/PHYLO/Materials/6_PHYLO_Seqtab jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiom__PHYLO/output


#   =====================================================================================================================================================


## 8?) UCHIME  ::  additional chimera removal
  # Usearch_64 as per SWALKER RY-suggestions. Do pre-RDSAgg!
  cd /data/jamie/PHYLO/Materials/6_PHYLO_Seqtab
  usearch_64 -uchime_ref PHYLO_run1_otus.fna -db /data/jamie/ref/chimera_slayer_gold_db/rRNA16S_gold_v20110519.fasta -uchimeout PHYLO_uc_out.uchime -strand plus -nonchimeras PHYLO_uc_good.fasta &> PHYLO_uc_usearch_chimera_removal.log


    ##
    ##  now need to extract those chimeras from ASV table!
    ##
    # all chimeric lines end with 'Y'
    grep 'Y$' PHYLO_uc_out.uchime | sed 's/\t9.*$//g' > ./PHYLO_uc_out_Y.uchime
    sed 's/\d*\.\d*\sSeq/faaa/' ./PHYLO_uc_out_Y.uchime > PHYLO_uchime_seqs.txt   # don't really need to do this
    # into R, and extact ASVs in that column...
     ~/R/bin/Rscript ~/PHYLO/analyses/8_UCHIME.stripChim.R      # work out docopt?




## 9?)  M O T H U R
  # question utility of this step, or at least provide a comparison. For now this (the proven approach) should remain the default.
  # OR FIND A NEW DEFAULT ( mockrobiome? )
  
  #fasta=/data/jamie/PHYLO/Materials/6_PHYLO_Seqtab/PHYLO_run1_otus_uc.fna,

  cd /data/jamie/PHYLO/Materials/6_PHYLO_Seqtab
  /home/jamie/BioApps/mothur/mothur '#classify.seqs(
  fasta=/data/jamie/PHYLO/Materials/6_PHYLO_Seqtab/PHYLO_run1_otus_uc.fna,
  reference=/home/jamie/BioApps/mothur/rdp/trainset16_022016.rdp.fasta,
  taxonomy=/home/jamie/BioApps/mothur/rdp/trainset16_022016.rdp.tax,
  processors=25,
  cutoff=80,
  probs=F)' &     # use fg to bring back to foreground
  
  sed 's/;/\t/g' PHYLO_run1_otus_uc.rdp.wang.taxonomy > PHYLO_MOTHUR_taxonomy.tsv


## 10?)   Q I I M E  :: phylogenetics using USEARCH output and Q1
    #   "There are a number of diversity metrics 
    #    like unifrac distance that require the 
    #    construction of a phylogenetic tree"  
    #              -  so true, so do.
    
    ## or do it in R using DECIPHER
    
  ## PARALLEL   :: update, parallelise ::  parallel_align_seqs_pynast.py -i PHYLO.run_otus.fna -t /usr/share/qiime/data/core_set_aligned.fasta.imputed -o ./done_parallel
  cd /data/jamie/PHYLO/Materials/6_PHYLO_Seqtab
  align_seqs.py -i PHYLO_run1_otus_uc.fna -t /usr/share/qiime/data/core_set_aligned.fasta.imputed
  filter_alignment.py -i ./pynast_aligned/PHYLO_run1_otus_uc_aligned.fasta
  make_phylogeny.py -i PHYLO_run1_otus_uc_aligned_pfiltered.fasta -r midpoint -o PHYLO_run1_uc_rooted.phylo.tre




# ==========================================================#     E   N   D   !    #============================================================



## Z) APPENDIX :: back up scripts as used.
  # breadcrumbs back in the box
  # tempting to rename & move output now - don't, in case not done correctly.   

  cp ./nohup.out ~/Materials
  cp ~/PHYLO/analyses/*.sh ~/Materials
  cp ~/PHYLO/analyses/*.R ~/Materials

  # scp -r /data/jamie/PHYLO/Materials/6_PHYLO_Seqtab jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiom__PHYLO/output/

# =====================================================#







##  =============   R E T H I N K   =============================================

## 5)  D A D A   2  ::   Deprecated -  Double check that this is a good idea.

  




##  =============   R E T I R E D   =============================================



## 0.2  ::   S A M P L E   Y O U R   R E A D S
  ## a silly thing: get a sample of reads to pass to learnErrors, should you not want to error-learn on all sample (or in WMI's case, have enormous samples that only need to be learned once.)
  # for i in $(ls *1.fastq.gz ) ; do zcat $i | sed -n '600001,1600000p;1600001q' >> /mnt/data/jfg_phd/seqweed/PEtrimmed_sample.R1.fastq ; done
  # for i in $(ls *2.fastq.gz ) ; do zcat $i | sed -n '600001,1600000p;1600001q' >> /mnt/data/jfg_phd/seqweed/PEtrimmed_sample.R2.fastq ; done
  # 
  # # recompress: gzip -c X > Y
  # gzip -c /mnt/data/jfg_phd/seqweed/PEtrimmed_sample.R1.fastq > /mnt/data/jfg_phd/seqweed/PEtrimmed_sample.R1.fastq.gz
  # gzip -c /mnt/data/jfg_phd/seqweed/PEtrimmed_sample.R2.fastq > /mnt/data/jfg_phd/seqweed/PEtrimmed_sample.R2.fastq.gz
  
### 0) different data commit scenarios
# scp ##################*gz /data/jamie/PHYLO/Materials/0_raw_reads   # on asp.ii
# ls ###############*gz | head -n 100 | xargs cp -t /data/jamie/PHYLO/Materials/0_raw_reads
# gunzip /data/jamie/PHYLO/Materials/0_raw_reads/*.gz     #not zipped
# ls jamie@143.239.189.162:/data/sidney/16s/biopsy/run1/*gz | head -n 100 | xargs cp -t /data/jamie/PHYLO/Materials/0_raw_reads
# ls /data/jamie/PHYLO/Materials/3_d2 | tail -n 500 | xargs mv -t /data/jamie/PHYLO/Materials/3_d2/3_D2_REST

  
### 4.4)  -  benched, as not enough reads used to estimate error model
# #    do echo ~/R/bin/Rscript ~/PHYLO/4.4_DADA2.error_pool.R --forwardReadpath d2.$i.R1.fastq --reverseReadpath d2.$i.R2.fastq --samplename $i.output --indir /data/jamie/PHYLO/Materials/3_d2 --outdir /data/jamie/PHYLO/Materials/4_error --samplelist ~/Materials/sample_run_list.txt >> ~/Materials/step_4_error.R  ;
# # for i in $(ls *.R1.fastq | cut -f 2 -d ".");     # catches sample ID, so can specify F/R
# #   do echo ~/R/bin/Rscript ~/PHYLO/analyses/4.4_DADA2.error_pool.R --forwardReadpath d2.$i.R1.fastq --reverseReadpath d2.$i.R2.fastq --samplename $i.output --indir /data/jamie/PHYLO/Materials/3_d2 --outdir /data/jamie/PHYLO/Materials/4_error --samplelist ~/Materials/sample_run_list.txt >> ~/Materials/step_4_error.R  ;
# # done 
# # 
# # parallel -j 16 < ~/Materials/step_4_error.R   2>&1 | tee -a ~/Materials/4_d2_error.log   # long hang at the end of this, presumably parallel shutting down
 

### 5) Average error for each run
# # based on a sample of samples-per-run (usually first 1x10^8 reads). Consider setting RANDOMIZE=TRUE for more representative -  are we still using this though? - trivial step to do, but a significant difference
# ## A S K
#   cd /data/jamie/PHYLO/Materials/4_error
#   ~/R/bin/Rscript ~/PHYLO/analyses/5.1_DADA2.avgError.R  2>&1 | tee -a ~/Materials/5_d2_avggErr.log


    

