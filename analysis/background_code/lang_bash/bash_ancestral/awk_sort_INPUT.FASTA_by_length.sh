#! /usr/bin/awk
# from https://www.biostars.org/p/153999/


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  INPUT.FASTA  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1n | cut -f 2- | tr "\t" "\n" |\
cat >> INPUT_SORTED.FASTA
