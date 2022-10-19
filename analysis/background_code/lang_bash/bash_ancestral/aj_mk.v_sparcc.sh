#! /usr/bin/zsh
cd /home/eri/jf/sparcc


## 1: Test OTU Table for Covariance, with default iterations and xiterations (10,20)
time python2.7 /home/eri/jf/sparcc/SparCC.py /home/eri/jf/sparcc/aj_sparcc/aj_otu.txt --cor_file=/home/eri/jf/sparcc/aj_sparcc/aj_cor.txt --cov_file=/home/eri/jf/sparcc/aj_sparcc/aj_cov.txt


## 2: Iterate a teech (1000) of shuffled OTU tables: (sorta fast!)
time python2.7 /home/eri/jf/sparcc/MakeBootstraps.py /home/eri/jf/sparcc/aj_sparcc/aj_otu.txt -n 1000 -t aj_permd_#.txt -p /home/eri/jf/sparcc/aj_sparcc/aj_boots/

cd aj_sparcc/aj_boots

## 3: SparCC on each bootstrap table as above:
time for i in $(ls *.txt); do python2.7 /home/eri/jf/sparcc/SparCC.py $i --cor_file cor_"$i" --cov_file cov_"$i"; done


## 4: Generate Pseudo-P Values for Correlation using Iterated, Correlated Datasets:
time python2.7 /home/eri/jf/sparcc/PseudoPvals.py '/home/eri/jf/sparcc/aj_sparcc/aj_cor.txt' '/home/eri/jf/sparcc/aj_sparcc/aj_boots/cor_aj_permd_1.txt' 1000 -o /home/eri/jf/sparcc/aj_sparcc/aj_one_sided_pp.txt -t one_sided 

