cd ~/jf/out/q2
mkdir ea
cd ea

##   T A X   T R A I N    F O R    E A
#already imported:
# s128_97_tax.qza

# create appropriate sequence training DB

	qiime feature-classifier extract-reads --i-sequences ../z.ref/s128_97_seq.qza --p-f-primer TAGATACCCSSGTAGTCC --p-r-primer CTGACGRCRGCCATGC --o-reads ../789-1068_ref-seqs.qza

# train

	qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ../789-1068_ref-seqs.qza --i-reference-taxonomy ../z.ref/s128_97_tax.qza --o-classifier ../789-1068_s128_classifier.qza




# CHECK MAPFILE
# less ea_map.txt

# FASTQ
gzip -k ../../0.start_out/ea.fastq ; mv ../../0.start_out/ea.fastq.gz sequences.fastq.gz 
gzip -k ../../1.splib_out_ea/barcodes.fastq ; mv ../../1.splib_out_ea/barcodes.fastq.gz barcodes.fastq.gz
mkdir ea_in ; mv *.gz ea_in


# IMPORT
qiime tools import --type EMPSingleEndSequences --input-path  ea_in/ --output-path ea_seqs.qza 


# DEMULTIPLEX
qiime demux emp-single --i-seqs ea_seqs.qza --m-barcodes-file ea_map.txt --m-barcodes-column BarcodeSequence --o-per-sample-sequences demux.qza
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
#qiime tools view demux.qzv

echo 'picking 300bp cutoff based on SILVA length summary'
# DADA2
qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left 20 --p-trunc-len 280 --p-n-threads 7 --o-representative-sequences rep-seqs.qza --o-table table.qza 


# TABLE
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file ea_map.txt
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
 

# TAXONOMY
qiime feature-classifier classify-sklearn --i-classifier ../789-1068_s128_classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file ea_map.txt --o-visualization taxa-bar-plots.qzv

# EXPORT
qiime tools export table.qza --output-dir exported-feature-table


