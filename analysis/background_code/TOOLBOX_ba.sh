# #---------------------------------------------------


##   I L L U M I N A   B A S E M O U N T  ##---------------------------------------------------
  
  ## see first the   C L I    O V E R V I E W 
  #
  #    https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
  # 

  ## an imaginary dir-type API for illumina data on AWS
    # needs to be configured (fill this bit in)
    
    # > basemount config ... ?      < ! ! ! > 
  
  ## and then mounted in the dir of your choice (one dir per user, all projects runs etc will be mirrored here)
    mkdir ~/data/ill_bsm
    basemount ~/data/ill_bsm
  
  ## can then ls, cp and etc from this point at leisure
  ## NOTES ::
    # initial access can be very slow, will accelerate as data provisioned (be patient!)
    # the api is natively parallel (def t=8), so don't recklessly stack parallel calls on it
  # e.g.
    ls -lsh ~/data/ill_bsm/Projects/SomeProj/Samples/*/Files/*gz                    # list all samples
    time cp ~/data/ill_bsm/Projects/SomeProj/Samples/samp_acc_ID/Files/*gz ~/data   # 13MB amplicons
      # > real	0m6.232s
      # > user	0m0.005s
      # > sys	0m0.027s




##   S U B S E T   F A S T Q S     ##---------------------------------------------------

  (
  OUT=/mnt/data/jamie/test
  SEQDIR=/mnt/data/jamie/metagenomes_017/
  cd $SEQDIR/
  # cd /mnt/workspace/jamie/dan/0_raw_test
  for i in $(ls *gz | head -10);
    do zcat $i | head -1000 > $OUT/test_${i/.gz/} ;
       gzip $OUT/test_${i/.gz/}
    done )
    


##  batches, bitches!               ##---------------------------------------------------
   (  BATCHES='2'
      FOUNT=$OLDH3
      QUANT=$(expr $(lk $FOUNT | grep -c 'gz') / $BATCHES )
        for i in $OLDH1 $OLDH3 ;  # $OLDH2   # just running 2x3j20t
          do 
            for j in $( ls $FOUNT | head -$QUANT );  do mv $FOUNT/$j ${i}/ ; done ;
          done
          )

  ## grab real sample names
      for i in $(cat $DAN/Materials/temp2); 
          do cd $KOUT ; ls *${i}* | sed 's/.kneaddata.fastq.gz//'  ; 
        done > $DAN/Materials/temp3



##   R E A D S   &   B A S E S            #---------------------------------------------------
  
    # https://bioinformatics.stackexchange.com/questions/935/fast-way-to-count-number-of-reads-and-number-of-bases-in-a-fastq-file
    fix_base_count() {
      local counts=($(cat))
      echo -e "${counts[0]}\t$((${counts[1]} - ${counts[0]}))"
    }
    export fix_base_counts

    fastq_summary(){
  
      ## assume passed a file name
        # sample=$(1/fastq.gz/)        # or     
        sample=$( ls $1 | sed -E 's/.*(Lib.*_[1,2]).*/\1/g' )
        
        size=$( du -sh $1 | sed -E 's/^([0-9]*.[0-9]*[A-Z]).*/\1/' )
          
        reads_bps=$( gzip -dc $1 \
                            | awk 'NR % 4 == 2' \
                            | wc -cl \
                            | fix_base_count )
                          
      echo -e "${sample}\t${size}\t${reads_bps}"
    }
    export fastq_summary
 
  # -------



##   G A R B L E D   F A S T Q    ##---------------------------------------------------

  # missing
    find /media/jfg/HDD_18-198/1_qc/*/*1_fastqc.zip | sed -E 's/.*(17[0-9]{4}).*/\1/' | sort > ~/m558_QCd_f
    find /media/jfg/HDD_18-198/1_qc/*/*2_fastqc.zip | sed -E 's/.*(17[0-9]{4}).*/\1/' | sort > ~/m558_QCd_r
    find /media/jfg/HDD_18-198/metagenomes_558/0_raw/*gz |  sed -E 's/.*(17[0-9]{4}).*/\1/' | sort | uniq > ~/m558_rFASTQ
    comm --total ~/m558_QCd_f ~/m558_rFASTQ
    comm --total ~/m558_QCd_r ~/m558_rFASTQ
  
    > ~/outstanding.txt
    for i in $(comm -13 ~/m558_QCd_f ~/m558_rFASTQ); do cd /media/jfg/HDD_18-198/metagenomes_558/0_raw/ ; ls *${i}*_1.fastq.gz >> ~/outstanding.txt ; done
    for i in $(comm -13 ~/m558_QCd_r ~/m558_rFASTQ); do cd /media/jfg/HDD_18-198/metagenomes_558/0_raw/ ; ls *${i}*_2.fastq.gz >> ~/outstanding.txt  ; done
    
    # double check no outputs exist  -  note different f/r FQ missing for different samples etc
    find /media/jfg/HDD_18-198/1_qc/* | grep '5156'

  # check damage  
  
    # eg Lib_Nextera_5156_172120_S21_2.fastq.gz
      less -S -p 5 /media/jfg/HDD_18-198/metagenomes_558/0_raw/Lib_Nextera_5156_172120_S21_2.fastq.gz   # not brooken?
      
    # faster?
    comm -13 ~/m558_QCd_f ~/m558_rFASTQ | 
        parallel -j 7 "if gzip -t /media/jfg/HDD_18-198/metagenomes_558/0_raw/*{}*_2.fastq.gz ; 
                              then
                                echo 'file is ok'
                              else 
                                echo 'file is corrupt'
                            fi"

    # check re-submitted files
    parallel -j 7 "if gzip -t {} ; 
                      then
                        echo 'file is ok'
                      else 
                        echo 'file is corrupt'
                    fi" ::: ./*fastq.gz
                    
    # all is good, PULL from vader
    parallel scp jfg@143.239.154.14:{} /mnt/data/jamie/metagenomes_017/ ::: $(ssh jfg@143.239.154.14 ls /home/jfg/Downloads/*fastq.gz)
    
                    


 ##   F A S T Q C   /   M U L T I Q C      #---------------------------------------------------
  
  # can define non-standard adapter sequences for FastQC, but gives a jlang excpetion and failure.
  # see ~/BioApps/fastqc/Configuration
  
    fastqc --threads 20 $KF/*fastq.gz -o .
    mkdir $DAN/1_qc/2_knead_f_qc ; mv $KF/*{zip,html} $DAN/1_qc/2_knead_f_qc
    
  # who knows
    conda install -n p27 -c bioconda multiqc
    multiqc $DAN/1_qc/2_knead_f_qc/
    
  # or wrapper  
    quick_qc(){
      local in=$1
      local out=$2
      local title=$3
      local threads=$4

      local work=$out/${title}_quick_qc
      
      mkdir $work
      cd $work
      fastqc -t $threads $in/*fastq.gz -o $work -f fastq  --noextract # like the custom adapter idea, but causing jlang exception atm
      
      multiqc -n $title $work
      rm $work/*fastqc.html
      cp $work/${title}.html $out
#      tar -czvf $out/${title}.tar.gz $work  -C $out $work   # cant work C and --remove-files deletes cwd 
      }
    export -f quick_qc




##   I D   F I L E S   O F   I N T E R E S T       #---------------------------------------------------

  get_id(){
    local dir=$1
    local term=*${2}*
      find $dir/$term | sed -r 's/.*[\/_]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{4,5}[AB]?$'
  }
  export get_id
  


##   Q   F I L E S   O F   I N T E R E S T       #---------------------------------------------------

  get_q(){
  # assume first folder=source, contains *gz files, and output is just outputs
    local inputs=$1
    local outputs=$2
      find $inputs/*gz | sed -r 's/.*[\/_]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{4,5}[AB]?$' > ./.get_q_temp_inputs
      find $outputs/*report | sed -r 's/.*[\/_]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{4,5}[AB]?$' > ./.get_q_temp_outputs
      comm -23 ./.get_q_temp_inputs ./.get_q_temp_outputs
      rm ./.get_q*
  }
  export get_q
  
  # should refactor so can pass any globbable expression
  
  ## eg, where both custom functions are declared 
  zcat $KF/Lib*$(get_q $KF $MOUT)*fastq.gz | less -S  
                    # | \
    # env_parallel meta25_banj $KF/Lib*_{}_*fastq.gz $MOUT




  
# ##   D O M E S T I C A T E          #---------------------------------------------------

  ## arguably v 0.4.0  -   domesticate version in shcrape 3.1~

# looks for matching h2.tar.gz (ref+srce), check if outputs exist, then H2's  

make_humann() {
    # better way of doing this surely
    if [[ $1 == "" ]]
    then
    echo "missing Args, quitting" ; exit 0
    fi
    
    # get_sample sample_name## $KF $HF/o_o $HOUT $SCRAPE 12
    
    local sample_name=$1      # make the name
    local input=$2
    local work=$3
    local outtar=$4
    local outtxt=$5
    local h2threads=$6
    
    local sample_path=$work/${sample_name}.kneaddata_humann2_temp
    
    # conda activate p27
    
    echo -e "\nID: $sample_name \nIN: $input \nWRK: $work \nTAR: $outtar \nTXT: $outtxt \nNAM: $sample_name \nPAT: $sample_path"
    
    # ====
    
    cd $work
    mkdir $sample_path
    
    if [ -f  $outtar/*$sample_name*tar.gz ]
    then 
    cp $( ls -t $outtar/*$sample_name*tar.gz | head -1 ) $work/   # copy first, newest match
      tar -xzvf $work/*$sample_name*.tar.gz --exclude='*/tmp*' --transform='s/.*\///' -C $sample_path/     # flatten output to -C
      fi
    
    
    if [ -f $sample_path/*$sample_name*abundance.tsv ] && [ -f $sample_path/*$sample_name*families.tsv ] && [ -f $sample_path/*$sample_name*coverage.tsv ]
    then 
    echo '+++ H U M A n N 2  :: f i l e s   e x i s t +++'
    cp $sample_path/$sample_name*{coverage,abundance,families,bugs}* $outtxt/
      rm -rf $work/*$sample_name*.tar.gz $work/*$sample_name*
      else
        echo -e "+++ no  H U M A n N 2  outputs - do humann2 on $sample_name +++\n+++\n+++\n+++ metaphlan2 DB should be chocophlan-v20 :: is this true? +++"
    humann2 --resume --threads $h2threads  --metaphlan-options='-x mpa_v20_m200 --bowtie2db /home/jamie/ref/chocophlan-v20' \
    --input $input/$sample_name.kneaddata.fastq.gz --memory-use maximum  --output $work &&
      cp $sample_path/*$sample_name*log* $outtxt/
      cp $work/*$sample_name*genefamilies* $outtxt/
      cp $work/*$sample_name*path* $outtxt/
      cp $sample_path/*$sample_name*bugs* $outtxt/
      scp $sample_path/*$sample_name*log* jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiome/SB__Dan/output/logs
    ssh jfg@143.239.154.14 echo "$sample_name done $(date)" >> /home/jfg/Dropbox/SeqBiome/SB__Dan/output/h2_humanised
    
    rm $work/*$sample_name*.tar.gz
    rm -rf $sample_path/tmp*
      tar -cv $work/*$sample_name* --remove-files | pigz -p $h2threads > $work/${sample_name}_hu2.dan2.tar.gz &&
      mv $work/${sample_name}_hu2.dan2.tar.gz $outtar/
      fi > $work/$sample_name.dan2.log 2>&1
    
}

export -f make_humann




  ## check multi logicals         #---------------------------------------------------
  
    if [ -f a.tsv ] && [ -f b.tsv ] && [ -f c.tsv ]
        then echo "files exist" ; echo $dir/*{a,b,c}*
    fi

  
  
  ## v0.1                         #---------------------------------------------------
  # get_sample() {
  #     local name=$( ls $2/*$1* | sed -r 's/.*(Lib.*_FR).*/\1/' )
  #     echo $name
  # }
  # export -f get_sample
  # 
  # get_sample 5492 $KF
  # parallel get_sample {} $KF :::  5492 5493 5670
  # parallel 'name=$(get_sample {} /mnt/data/jamie/dan_fin/2_knead_f) ; echo $name' :::  5492 5493 5670
  # lk $HF/outstanding_outlier/*_hu2.tar.gz | sed -r 's/.*[\/_]([0-9]{4}[A-Z]?)[\/_].*hu2.tar.gz/\1/' | tail -3 |  \
  #         parallel -j 3 sample_name=$( get_sample {} $KF ) \; $sample_name 


# --------------------------------------------------------------------------------------------------------------
    
    feed_knead(){
    
      # better way of doing this surely
      if [[ $1 == "" ]]
      then
      echo "missing Args, quitting" ; exit 0
      fi
    
      # ---    
    
     local seq=$1
     local kin=$2
     local kout=$3
     local db=$4
     local procs=$5
      
      # check
      echo -e "__using variables:_______________\
      \nInput sample ID :: $seq \nKnead input :: $kin \
      \nKnead output :: $kout \nKneaddata Database :: $db\
      \nthreads per knead job :: $procs \n________________________________"
      
      # ---
      
      # KNEADDATA
      time kneaddata -t ${procs} -i $kin/${seq}_1.fastq.gz -i $kin/${seq}_2.fastq.gz \
      -o $kout/$seq -db $db --max-memory 45g --remove-intermediate-output \
      --trimmomatic ~/.conda/envs/p27/share/trimmomatic --trimmomatic-options="ILLUMINACLIP:NexteraPE-PE.fa:2:30:10" \
      --trimmomatic-options="CROP:145" --trimmomatic-options="HEADCROP:15" --trimmomatic-options="SLIDINGWINDOW:4:20" \
      --trimmomatic-options="MINLEN:50" > $kout/knead_${seq}.stout 2>&1
      
      ##  overheads
      seqfr=${seq}_FR.kneaddata.fastq                 ## paired end
      
      ## stow quarantined reads
      mkdir $kout/$seq/quarantine
      mv $kout/$seq/*contam* $kout/$seq/quarantine
      tar -czvf $kout/$seq/${seq}_quarantine.tar.gz $kout/$seq/quarantine --remove-files &
      
      ## cat, gz & then delete usable reads
      mkdir $kout/$seq/good_reads
      mv $kout/$seq/*fastq $kout/$seq/good_reads
      cat $kout/$seq/good_reads/* > $kout/$seq/$seqfr
      gzip $kout/$seq/$seqfr
      
      rm -rf $kout/$seq/good_reads                            # deletes intermediate fastqs
    
      }
      
      export -f feed_knead




  
##   T A X   C O M P A R E S         ##---------------------------------------------------
    
  # from kraken
  
  gettax_krak(){
    
      local code=$1
      local ref=$2
      
        if [ $code -eq "0" ]
        then
          echo -e "0\t|\tundefined\t|\t|\t unname name |"
        else
          awk -F '\t' '$1 == "'"$code"'"' $ref | grep 'scientific name'          # !!! note awk!
        fi
      }
    export -f gettax_krak 

  # then, but super slow:    
    for i in $( head bracken_444_combo | cut -f 2 ); do gettax_krak $i $DB/kraken2_stnd/taxonomy/names.dmp ; done

  
  ## from motus
      
  get_tax_motus(){
    
    local metam=$1
    
    mgc=$(grep $metam ~/BioApps/m2/db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv | cut -f 4 )
    motu=$(grep $mgc ~/BioApps/m2/db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv | cut -f 1  | grep 'mOTU' )
    tax=$( grep $motu ~/BioApps/m2/db_mOTU/db_mOTU_taxonomy_CAMI.tsv | sed -E 's/(^.*_mOTU_v25_[0-9]{5}).*\|(.*)$/\1 \2/')
    
    echo -e "$1\t$mgc\t$motu\t$tax"
    }
  export -f get_tax_motus


  ## get uniq fastq headers from file 
  
    u_head_out(){
      local file=$1
      # echo "get fastq headers from $file with str A00.*\#"     # gets printed to file! XD
        sed -E 's/(A00.*\#0\/[1,2]).*/\1/' $file | sort | uniq
        # sed -E 's/(A00.*\#0).*/\1/' $file | sort | uniq
      }
    export -f u_head_out


  # -------

  ## is integer?
  for i in $()
    [[ $a =~ ^-?[0-9]+$ ]] && echo integer


  

# --------------------------------------------------------------------------------------------------------------

  ## better metaphlan - do v.25 with sensible F/R 
  # currently banjaxed
  
meta25_banj(){

  local knead=$1
  local mout=$2
  local met_name=$(basename $knead)
  local met_name=${met_name/_FR.kneaddata.fastq.gz/}
  
      ## re-partition F&R reads
        # local split_f=${met_name}_F.fastq
        # local split_r=${met_name}_R.fastq
      ## problem reading split fastqs
        # zgrep -A 3 '#0\/1' $knead > $mout/$split_f  &&
        # zgrep -A 3 '#0\/2' $knead > $mout/$split_r  &&
      ## problem here could be the renaming, but prob not - the following neither works nor makes a difference
        # parallel "name=$(ls  ./{} | sed 's/_FR.kneaddata.fastq.gz//' ) ; zgrep -A 3 '#0\/1' {} > ${name}_F.kd.fastq" ::: ./*fastq.gz &
        # rm $mout/$split_f &
        # rm $mout/$split_r &

    metaphlan2.py $knead --input_type fastq -t rel_ab_w_read_stats --bowtie2out $mout/${met_name}.bt2 -o $mout/${met_name}_metap25 --nproc 10 --force
}
export meta25_banj
  
  env_parallel -j 7 meta25 {} $MOUT ::: $KF/*gz
  
  

# ---------------------------------------------------------------------------------------------------------

      # first awk!
        # awk -F '\t' '$3 == "d"' filename
        # awk '/search_pattern/ { action_to_take_on_matches; another_action; }' file_to_parse

      # from a file with codes, lookup in the k2 db
        gettax(){
          local code=$1
          local ref=$2
          
            if [ $code -eq "0" ]
            then
              echo -e "0\t|\tundefined\t|\t|\t unname name |"
            else
              # https://www.biostars.org/p/97653/    thanks Ole
              awk -F '\t' '$1 == "'"$code"'"' $ref | grep 'scientific name'
            fi
          }
        export -f gettax 

        parallel gettax {} $DB/kraken2_stnd/taxonomy/names.dmp :::: $sub/codes.tmp  > $sub/names.txt


# ---------------------------------------------------------------------------------------------------------




## bash hints

	# preserve your params - curly brackets: {}
		${i}_wont_contage


	# parameter modification
		${parameter/pattern/string} — replace the first match of pattern with string.
		${parameter//pattern/string} — replace all matches of pattern with string.

	# call multiples in parallel
	parallel -j 2 test.sh {1} {2} $KOUT $HOUT ::: $READS/*.1.fastq ::: $READS/*.2.fastq

	## multicrush in parallel (all cores)
	  # note --remove-files is at the end  -  f flag expects next arg to be the name of archive created
		parallel tar -czvf {}.tar.gz {} --remove-files ::: ./*  


	# sticking assigned vars inside scripts wont work if you send them to someone else to process! 
	# esp true if that someone is parallel etc. 
	#i.e.
		cat sniff.txt > ./knead.log 
	# not 
		cat stiff.txt > $MISTER/knead.log


	## multiple lines
	  this.sh --will work /
	    --just-fine
	  # as there's _nothing_ after the backslash
	  but_this.sh --will work /  # not at all
	    --because-of-the-#d-portion --not-being-an-arg   # so watch out
  







