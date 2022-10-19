## KO_reformat :: jfg v0.2 Aug'19

# from https://www.genome.jp/kegg/ko.html, go "KEGG Orthology (KO)", i.e. https://www.genome.jp/kegg-bin/get_htext?ko00001, then "Download htext"
# file saved and modified (HTML header/footer removed) from https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir= ,  2019.07.03

k <- data.table::fread('~/Dropbox/SeqBiome/sb_4.11/analysis/background_code/ko_jfg_tsv.keg', sep='\t', header=FALSE, fill = TRUE)    # problems.com with characters and reading in

A <- 'bumfluff'
for (i in 1:nrow(k)){ if (k[i,2] ==''){k[i,2] <- A}else{A <- k[i,2]} }         # slooow loops, populate indents
for (i in 1:nrow(k)){ if (k[i,3] ==''){k[i,3] <- A}else{A <- k[i,3]} }         # apply(k2, 2, function(x){ if (x ==''){x <- A}else{A <- x} })
for (i in 1:nrow(k)){ if (k[i,4] ==''){k[i,4] <- A}else{A <- k[i,4]} }

k <- k[k$V1 == 'D' ,]
colnames(k) <- c('let', 'A', 'B', 'C', 'D', 'EC')
knitr::kable(head(k))

# when you're doing things once only ever, programmatic is a luxury timespend
A <- stringr::str_split_fixed(k$A , n=2, ' ')
B <- stringr::str_split_fixed(k$B , n=2, ' ')
C <- stringr::str_split_fixed(k$C , n=2, ' ') ; C1 <- C[,1] ; C2 <- cbind(stringr::str_split_fixed(C[,2] , n=2, ' \\[PATH:'))  
D <- stringr::str_split_fixed(k$D , n=2, ' ')
EC <- k[,'EC']

k <- as.matrix(cbind(A,B,C1,C2,D, EC))
k <- gsub('\\[|\\]', '', k)    # this deletes your computer
colnames(k) <- c('A_num' , 'A_categ', 'B_num' , 'B_categ', 'C_num' , 'C_categ', 'D_PATH' , 'D_ORTH', 'D_product' , 'D_EC' )
rownames(k) <- paste(k[, 'D_PATH'],k[, 'D_ORTH'], 1:nrow(k), sep='_')
# knitr::kable(k[1:10,])

saveRDS(k, '~/Dropbox/SeqBiome/sb_4.11/analysis/background_code/r__KO_hierarchy.RDS')
write.table(k, '~/Dropbox/SeqBiome/sb_4.11/analysis/background_code/ko_jfg_fullformat.tsv', sep='\t')

## returns DF of the form:   -  see gene ID's in second-last column that might by of use?

    # A_num   A_categ      B_num   B_categ                   C_num   C_categ                        D_PATH    D_ORTH   D_product                                                 D_EC                
    # ko00010_K00844_1 "09100" "Metabolism" "09101" "Carbohydrate metabolism" "00010" "Glycolysis / Gluconeogenesis" "ko00010" "K00844" " HK; hexokinase"                                         "EC:2.7.1.1"        
    # ko00010_K12407_2 "09100" "Metabolism" "09101" "Carbohydrate metabolism" "00010" "Glycolysis / Gluconeogenesis" "ko00010" "K12407" " GCK; glucokinase"                                       "EC:2.7.1.2"        
    # ko00010_K00845_3 "09100" "Metabolism" "09101" "Carbohydrate metabolism" "00010" "Glycolysis / Gluconeogenesis" "ko00010" "K00845" " glk; glucokinase"                                       "EC:2.7.1.2"        
    # ko00010_K01810_4 "09100" "Metabolism" "09101" "Carbohydrate metabolism" "00010" "Glycolysis / Gluconeogenesis" "ko00010" "K01810" " GPI, pgi; glucose-6-phosphate isomerase"                "EC:5.3.1.9"        
    # ko00010_K06859_5 "09100" "Metabolism" "09101" "Carbohydrate metabolism" "00010" "Glycolysis / Gluconeogenesis" "ko00010" "K06859" " pgi1; glucose-6-phosphate isomerase, archaeal"          "EC:5.3.1.9"        
    # ko00010_K13810_6 "09100" "Metabolism" "09101" "Carbohydrate metabolism" "00010" "Glycolysis / Gluconeogenesis" "ko00010" "K13810" " tal-pgi; transaldolase / glucose-6-phosphate isomerase" "EC:2.2.1.2 5.3.1.9"
