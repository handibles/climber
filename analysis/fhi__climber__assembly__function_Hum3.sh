

##   SB__Dan_extreme_catchup

## Infrastruc

    DSK=jfg@143.239.154.14:/home/jfg/Dropbox/SeqBiome/SB__Dan/output/
      
    DATA=/mnt/data/jamie                  # VADER
    RAW=/mnt/data/jamie/dan/0_raw         # VADER
    DB=/home/jamie/ref                    # VADER
    DAN=/mnt/workspace/jamie/dan3          # VADER         !!! new DAN dir  !!!

    READS=$DAN/0_raw
    QC=$DAN/1_qc
    KOUT=$DAN/2_kneaded
    MOUT=$DAN/3_m2
    HOUT=$DAN/4_h2
    SCRAPE=$DAN/5_scrape
    
    # UR90DB=$DB/uniref
    conda activate p27    
    
    # check here that all true
    mkdir $DAN
    mkdir $DAN/Materials
    mkdir $READS
    mkdir $QC
    mkdir $KOUT
    mkdir $MOUT
    mkdir $HOUT
    mkdir $SCRAPE
    
    KF=/mnt/data/jamie/dan_fin/2_knead_f
    MF=/mnt/data/jamie/dan_fin/3_m2_f
    HF=/mnt/data/jamie/dan_fin/4_h2_f
    SF=/mnt/data/jamie/dan_fin/5_scrape_f

  
  
  
# ===== #   F U N C T I O N S  # ================================================ #

  ## source the FNs   
    # # or...
    #   scp ~/Dropbox/SeqBiome/sb_4.11_master/analysis/background_code/TOOLBOX_ba.sh jamie@143.239.186.42:/home/jamie/
    #   . ~/TOOLBOX_ba.sh
  
  
    
# ================================================= #  --------------------------------------------------------------------------------------------------------------------
# ==== #   G E T   T  O   I T  # ================== #  --------------------------------------------------------------------------------------------------------------------
# ================================================= #  --------------------------------------------------------------------------------------------------------------------

## Checking

    ## dead, dont worry: 172588, 172269, 172270

    ## need to mend some samples through pipeline
    # 172123
    # 172125
    # 172127
    # 172128
    # 172129
    # 172130
    # 172327
    # 172328

  # sample_id
    # 4-char code
    # aws s3 ls s3://danhout | sed -r 's/.*[\/_ ]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq > $DAN/Materials/tempA    # extra space in [there]
    find $DATA/metagenomes_558/*gz | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp0
    find $KF/*gz | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp2
    find $MF/ | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp3
    find $HF/*gz | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp4
    find $HF | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp4all
    find $HF/o* | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp4.o
    find $SF/p.cov | sed -r 's/.*[\/_]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq  | grep -E '^[0-9]{4,5}[AB]?$' > $DAN/Materials/temp5.c
    find $SF/p.ab | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$'  > $DAN/Materials/temp5.a
    find $SF/gene | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp5.g
    find $SF/bugs/*txt.gz | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp5profs
    find $SF/bugs/*tsv.gz | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp5bugs
    find $KRK/*bracken | sed -r 's/.*[\/_]([0-9]{6})[\/_].*/\1/g' | sort | uniq | grep -E '^[0-9]{6}$' > $DAN/Materials/temp_br
    
    wc -l $DAN/Materials/temp*
    comm -23 $DAN/Materials/temp0 $DAN/Materials/temp2 

  # should have a step to check have all the files  -   as you didn't


## check missing files
    comm -23 $DAN/Materials/temp0 $DAN/Materials/temp2 | \
    parallel -j 8 "if gzip -t $DATA/metagenomes_558/*{}*gz ; 
                      then
                        echo 'file is ok'
                      else 
                        echo 'file is corrupt'
                    fi" 
    
    ## all is good, but we're staying on this server..
    # parallel rsync -avh jfg@143.239.154.14:{} /mnt/data/jamie/metagenomes_017/ ::: $(ssh jfg@143.239.154.14 ls /home/jfg/Downloads/*fastq.gz)
    
    # 008 just cp from source   -    017 files moved differently
    comm -23 $DAN/Materials/temp0 $DAN/Materials/temp2 | parallel cp $DATA/metagenomes_558/*{}*  $READS/ 

  
## Kneading ## ------------------------------------------------------
  
  # define 
    kinput=$READS
    koutput=$KOUT
    db=$DB/kneads
    kthreads=5
    
  # run

    ls $kinput/*fastq.gz  | sed -E 's/.*(Lib.*_S[0-9]{0,3}).*/\1/g' | sort | uniq |  \
        parallel -j 10 feed_knead {} $kinput  $koutput $DB/kneads 5

  # mv all fastqs into the open  
    mv $KOUT/*/*kneaddata.fastq.gz $KOUT
    lk $KOUT    
    
    
##   Re HUMANNising ## -----------------------------------------------

     #  local sample_name=$1      # make the name
     #  local input=$2
     #  local work=$3
     #  local outtar=$4
     #  local outtxt=$5
     #  local h2threads=$6

    # define 
      # {}={}           # 1
      hinput=$KOUT      # 2  # ref
      hwork=$HOUT       # 3
      tar_dir=$HOUT     # 4 # tar srce/dest
      houttxt=$SCRAPE   # 5  # for h2 gene/path files
      h2threads=15      # 6

      # db=$DB    # set as chocophlan-v20

    ## jamie....  
      # find $KOUT | grep 'kneaddata.fastq.gz' | sed -r 's/.*[\/_]([0-9]{4}[AB]?)[\/_].*/\1/g' | sort | uniq | head- 1
      #   parallel -j 4 make_humann {} $hinput $hwork $tar_dir $houttxt $hthreads #>> $HOUT/h2_domestest_basic.log 2>&1 

    # escape quotes!
      parallel -j 4 humann2 --threads $h2threads --metaphlan-options=\"-x mpa_v20_m200 --bowtie2db /home/jamie/ref/chocophlan-v20\" \
          --input {} --memory-use maximum  --output $hwork > $HOUT/{}_extra.log ::: \
          $KOUT/*fastq.gz >> $HOUT/rederehumann.log
    
    
    