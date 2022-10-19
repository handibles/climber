cd ~/jf/out/q2
mkdir dw
cd dw

# less dw_map.txt

# FASTQ
gzip -k ../../0.start_out/dw.fastq ; mv ../../0.start_out/dw.fastq.gz sequences.fastq.gz 
gzip -k ../../1.splib_out_dw/barcodes.fastq ; mv ../../1.splib_out_dw/barcodes.fastq.gz barcodes.fastq.gz
mkdir dw_in ; mv *.gz dw_in


# IMPORT
qiime tools import --type EMPSingleEndSequences --input-path  dw_in/ --output-path dw_seqs.qza 


# DEMULTIPLEX
qiime demux emp-single --i-seqs dw_seqs.qza --m-barcodes-file dw_map.txt --m-barcodes-column BarcodeSequence --o-per-sample-sequences demux.qza
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
#qiime tools view demux.qzv
echo 'picking cutoff@355bp based on Demult viz and shitty results @500'

# DADA2 - trunc@355
qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left 20 --p-trunc-len 355 --p-n-threads 7 --o-representative-sequences rep-seqs.qza --o-table table.qza


# TABLE
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file dw_map.txt
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
 

# TAXONOMY
qiime feature-classifier classify-sklearn --i-classifier ../905-1492_s128_classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file dw_map.txt --o-visualization taxa-bar-plots.qzv


# EXPORT
qiime tools export table.qza --output-dir exported-feature-table


