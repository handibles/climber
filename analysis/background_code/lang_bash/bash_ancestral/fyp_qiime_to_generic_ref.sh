
# ## from 
# >SH497095.07FU_LC146734_reps_singleton
# TTGAATTTCACCGGTTTGGTCTGTTGCTGGTTCCGAAAGGTTCATGTGCACGCCTCGCCTCTGATATCTCACCACCTGTGAACTTTAGTGGGCTGTGACGGCCTTGTCTTCAGCAGTTTCGTGTTGGGTTTTGGGTACTTGTACCTGCCAATGCGCTCTGCTGGG...
# ## and 
# SH497095.07FU_LC146734_reps_singleton   k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Tricholomataceae;g__Flagelloscypha;s__Flagelloscypha_japonica

# ## to
# >Boletales_sp|UDB004660|SH465869.07FU|reps_singleton|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Boletales;f__unidentified;g__unidentified;s__Boletales_sp
# GGAAGGACATTATCGAAACAAATGGGGGGAAGACTGTCGCTGGCCCTCGGGCATGTGCACGTCGACCTCTTCATACACACACACCTGTGCACCTTTGGTAGGTCTTCGAAAGAGGATCTATGATTATCATCACACCCTGTCGTATGGCCAGAATGTCTATATCAC...


# =====  M K 1 . 0 =========#
# 
## preferred method is PROCESS SUBSITUTION, as does not muck with env.vars
# https://wiki.bash-hackers.org/syntax/expansion/proc_subst

while read i;
  do echo '>'$i | sed 's/\sk__/|k__/g' >> new_ref.txt ; 
  echo $i | cut -f 1 -d \s | grep --file="-" -A 1 min_97.fasta | sed -n 2p >> new_ref.txt ;
done <  min_97.txt

  # ## gives:
  # >SH497095.07FU_LC146734_reps_singleton|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Tricholomataceae;g__Flagelloscypha;s__Flagelloscypha_japonica
  # TTGAATTTCACCGGTTTGGTCTGTTGCTGGTTCCGAAAGGTTCATGTGCACGCCTCGCCTCTGATATCTCACCACCTGTGAACTTTAGTGGGCTGTGACGGCCTTGTCTTCAGCAGTTTCGTGTTGGGTTTTGGGTACTTGTACCTGCCAATGCGCTCTGCTGGGGGG...


## =====  M K 1 . 2 =========#
  # need:
  # >Boletales_sp|UDB004660|SH465869.07FU|reps_singleton|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Boletales;f__unidentified;g__unidentified;s__Boletales_sp
  # >Sp_handle | ACC1 | ACC2 | status | taxonomy
  # [tax field 2.1:6]  [ acc2 ] [ acc 1 ] [ ]

while read i;
  do echo '>'$i | sed 's/\sk__/|k__/g' >> new_ref.txt ; 
  echo $i | cut -f 1 -d \s | grep --file="-" -A 1 97.fasta | sed -n 2p >> new_ref.txt ;
done <  97.txt



# ============

## fucking for-loops man.
# 
# IFS=$'\n'                       # set env value, as for loop will split on any whitespace
# 
# rm new_ref.txt
# for i in $(cat min_97.txt); 
#   do echo '>'$i | sed 's/\sk__/|k__/g' >> new_ref.txt ; 
#   echo $i | cut -f 1 -d \s | grep --file="-" -A 1 min_97.fasta | sed -n 2p >> new_ref.txt ;
# done 
# 
# unset IFS                       # fix env value


