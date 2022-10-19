## selfmade


##   F A S T A    T O    F A S T Q ,   W . /   B A R C O D E S 
convert_fastaqual_fastq.py -f 0.start_out/aj.fasta -q 0.start_out/aj.qual -F -b -o 0.start_out/ -c fastaqual_to_fastq
# g e t   b a r c o d e s
extract_barcodes.py -f 0.start_out/aj.fastq -o 1.splib_out_aj/ -c barcode_single_end -l 10  # woo!
# see r script for Splib Details


	# EA
	convert_fastaqual_fastq.py -f 0.start_out/ea.fasta -q 0.start_out/ea.qual -F -b -o 0.start_out/ -c fastaqual_to_fastq
	extract_barcodes.py -f 0.start_out/ea.fastq -o 1.splib_out_ea/ -c barcode_single_end -l 10  

	# DW
	convert_fastaqual_fastq.py -f 0.start_out/dw.fasta -q 0.start_out/dw.qual -F -b -o 0.start_out/ -c fastaqual_to_fastq
	extract_barcodes.py -f 0.start_out/dw.fastq -o 1.splib_out_dw/ -c barcode_single_end -l 10  




##   T A X   T R A I N
##
# tax training

	wget -O "s128.tgz" "https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_release.tgz"
	tar -xvzf s128.tgz
		
	qiime tools import --type 'FeatureData[Sequence]' --input-path s128_97.fasta --output-path s128_97_seq.qza
	qiime tools import --type 'FeatureData[Taxonomy]' --source-format HeaderlessTSVTaxonomyFormat --input-path s128_97.txt --output-path s128_97_tax.qza

# create appropriate sequence training DB

	qiime feature-classifier extract-reads --i-sequences s128_97_seq.qza --p-f-primer TGAAACTYAAAGGAATTG --p-r-primer GGTTACCTTGTTACGACTT --o-reads dwaj/905-1492_ref-seqs.qza

# train

	qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads dwaj/905-1492_ref-seqs.qza --i-reference-taxonomy s128_97_tax.qza --o-classifier dwaj/905-1492_s128_classifier.qza


			## havent created an appropriate test dataset, so mothball
			## test on eg dataset
			#qiime feature-classifier classify-sklearn \
			#  --i-classifier dwaj/905-1492_s128_classifier.qza \
			#  --i-reads dwaj/dwaj_rep-seqs.qza \
			#  --o-classification dwaj/dwaj_taxonomy.qza
			#qiime metadata tabulate \
			#  --m-input-file dwaj/dwaj_taxonomy.qza \
			#  --o-visualization dwaj/dwaj_taxonomy.qzv


#   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   
#==============================================================================================
#   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   =   


##
##   Big fasta of seqs and a big fastq of barcodes  =  E M P r o t o c o l
##
##   https://docs.qiime2.org/2018.4/tutorials/moving-pictures/


## keep it simple, keep it safe

	mkdir aj
	cd aj


# check map   (how did that get here if we just created the folder?...)

	less aj_map.txt


# need to imprt DW, AJ, and EA separately however (barcode conflict). Can do this in three diff setions, as not parallel'd, and 
# barcodes and filtering are carried out separately
	
	gzip -k ../../0.start_out/aj.fastq ; mv ../../0.start_out/aj.fastq.gz sequences.fastq.gz 
	gzip -k ../../1.splib_out_aj/barcodes.fastq ; mv ../../1.splib_out_aj/barcodes.fastq.gz barcodes.fastq.gz
	mkdir aj_in ; mv *.gz aj_in


## IMPORT 
# (EMProtocol, i.e. dir with barcodes and seqs in fastq.gz)

	qiime tools import --type EMPSingleEndSequences --input-path  aj_in/ --output-path aj_seqs.qza 


## DEMULT:

	qiime demux emp-single --i-seqs aj_seqs.qza --m-barcodes-file aj_map.txt --m-barcodes-column BarcodeSequence --o-per-sample-sequences demux.qza
	qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
	qiime tools view demux.qzv
	# trim first 30 bp to 'definitely' exclude primers etc?


## DENOISE
# Len trunc is critical: longer will be cut, shorter will be discarded: all seqs will be this length
# based on demux plot & silva, try 355 / 280 
	qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left 20 --p-trunc-len 355 --p-n-threads 7 --o-representative-sequences rep-seqs.qza --o-table table.qza

	#  --o-denoising-stats denoise_stats.qza  #deprecated?


## TABLE

	qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file aj_map.txt
	qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
 

## TAXONOMY
# best classify your taxonomy assignment on yourdwaj data...
# https://docs.qiime2.org/2018.4/tutorials/feature-classifier/

	qiime feature-classifier classify-sklearn --i-classifier ../905-1492_s128_classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
	qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

	qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file aj_map.txt --o-visualization taxa-bar-plots.qzv

## EXPORT

	qiime tools export table.qza --output-dir exported-feature-table



##   T R E E 
# phylo tree for other uses (UNIFRAC): Align, mask, tree, mid-root

	qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza
  	qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
	qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
	qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

 
##   A L P H A / B E T A   D I V 
## diversity metrics
# sampling depths based on group of second-sparsest samples, dropping 'sig' lowest  

	 qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 1109 --m-metadata-file  aj_map.txt --output-dir core-metrics-results

# Alpha
	qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file  aj_map.txt --o-visualization core-metrics-results/faith-pd-group-significance.qzv
	qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file  aj_map.txt --o-visualization core-metrics-results/evenness-group-significance.qzv

# Beta
	qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file  aj_map.txt --m-metadata-column Reactor --o-visualization core-metrics-results/unweighted-unifrac-reactor-significance.qzv --p-pairwise
  
# diverstiy and distance: hypothesis factory
	qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file  aj_map.txt  --o-visualization core-metrics-results/uw-unifrac-emperor-reactor.qzv
#--p-custom-axes Reactor

#Alph-rare
	qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 4000 --m-metadata-file  aj_map.txt --o-visualization alpha-rarefaction.qzv


