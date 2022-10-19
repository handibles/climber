#!/bin/bash

# n'oublie pas: in the sample data, OTUs are the rows. 

cd /home/jamie/BioApps/sparcc

	## 1: Test OTU Table for Covariance, with default iterations and xiterations (10,20)
	time python2.7 /home/jamie/BioApps/sparcc/SparCC.py /home/jamie/BioApps/sparcc/ry_sparcc/ry100_otu.txt --cor_file=/home/jamie/BioApps/sparcc/ry_sparcc/ry_cor.txt --cov_file=/home/jamie/BioApps/sparcc/ry_sparcc/ry_cov.txt

	## 2: Iterate a teech (15) of shuffled OTU tables: (fast!)
	time python2.7 /home/jamie/BioApps/sparcc/MakeBootstraps.py /home/jamie/BioApps/sparcc/ry_sparcc/ry100_otu.txt -n 500 -t ry_permd_#.txt -p /home/jamie/BioApps/sparcc/ry_sparcc/ry_boots/

cd ry_sparcc/ry_boots

	## 3: SparCC on each bootstrap table as above:
	## this script making is ugly, gotta be a better way
	for i in $(ls *.txt); 
		do echo python2.7 /home/jamie/BioApps/sparcc/SparCC.py $i --cor_file cor_"$i" --cov_file cov_"$i" ;
	done > ry_bootstrap.sh
	time parallel -j 26 < ./ry_bootstrap.sh
 

	## W O R K S   I F   done within py2.7 venv on bio-linux 
	## 4: Generate Pseudo-P Values for Correlation using Iterated, Correlated Datasets:
	time python2.7 /home/jamie/BioApps/sparcc/PseudoPvals.py '/home/jamie/BioApps/sparcc/ry_sparcc/ry_cor.txt' '/home/jamie/BioApps/sparcc/ry_sparcc/ry_boots/cor_ry_permd_#.txt' 500 -o /home/jamie/BioApps/sparcc/ry_sparcc/ry_one_sided_pp.txt -t one_sided 

	
	
	## test accession not passed to iterated tables
	# 2: 	# do for pearson also
	cd /home/jamie/BioApps/sparcc/ry_sparcc
	mkdir /home/jamie/BioApps/sparcc/ry_sparcc/ry_pear_boots
	python2.7 /home/jamie/BioApps/sparcc/SparCC.py /home/jamie/BioApps/sparcc/ry_sparcc/ry100_otu.txt --cor_file=/home/jamie/BioApps/sparcc/ry_sparcc/ry_pear_cor.txt --cov_file=/home/jamie/BioApps/sparcc/ry_sparcc/ry_pear_cov.txt -a pearson
	python2.7 /home/jamie/BioApps/sparcc/MakeBootstraps.py /home/jamie/BioApps/sparcc/ry_sparcc/ry100_otu.txt -n 500 -t ry_permd_#.txt -p /home/jamie/BioApps/sparcc/ry_sparcc/ry_pear_boots/
	cd /home/jamie/BioApps/sparcc/ry_sparcc/ry_pear_boots
## problems here -windows related?
	for i in $(ls *.txt) ; do echo python2.7 /home/jamie/BioApps/sparcc/SparCC.py $i --cor_file cor_"$i" --cov_file cov_"$i" ; done > ry_pear_bootstrap.sh
	parallel -j 26 < ./ry_pear_bootstrap.sh
## source venv here
source ~/BioApps/ry_venv/bin/activate
	time python2.7 /home/jamie/BioApps/sparcc/PseudoPvals.py '/home/jamie/BioApps/sparcc/ry_sparcc/ry_pear_cor.txt' '/home/jamie/BioApps/sparcc/ry_sparcc/ry_pear_boots/cor_ry_permd_#.txt' 500 -o /home/jamie/BioApps/sparcc/ry_sparcc/ry_pear_one_sided_pp.txt -t one_sided 
# have to unset here for the rest of the work
deactivate
	
	
	#3: 	# aaand spearman
	cd /home/jamie/BioApps/sparcc/ry_sparcc
	mkdir /home/jamie/BioApps/sparcc/ry_sparcc/ry_spe_boots
	python2.7 /home/jamie/BioApps/sparcc/SparCC.py /home/jamie/BioApps/sparcc/ry_sparcc/ry100_otu.txt --cor_file=/home/jamie/BioApps/sparcc/ry_sparcc/ry_spe_cor.txt --cov_file=/home/jamie/BioApps/sparcc/ry_sparcc/ry_spe_cov.txt -a spearman
	python2.7 /home/jamie/BioApps/sparcc/MakeBootstraps.py /home/jamie/BioApps/sparcc/ry_sparcc/ry100_otu.txt -n 500 -t ry_permd_#.txt -p /home/jamie/BioApps/sparcc/ry_sparcc/ry_spe_boots/
## problems here -windows related?
	cd /home/jamie/BioApps/sparcc/ry_sparcc/ry_spe_boots
	for i in $(ls *.txt); 
		do echo python2.7 /home/jamie/BioApps/sparcc/SparCC.py $i --cor_file cor_"$i" --cov_file cov_"$i" ;
	done > ry_spe_bootstrap.sh
	parallel -j 26 < ./ry_spe_bootstrap.sh
## source venv here
source ~/BioApps/ry_venv/bin/activate
	python2.7 /home/jamie/BioApps/sparcc/PseudoPvals.py '/home/jamie/BioApps/sparcc/ry_sparcc/ry_spe_cor.txt' '/home/jamie/BioApps/sparcc/ry_sparcc/ry_spe_boots/cor_ry_permd_#.txt' 500 -o /home/jamie/BioApps/sparcc/ry_sparcc/ry_spe_one_sided_pp.txt -t one_sided 
deactivate
