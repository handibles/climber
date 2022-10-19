##   N O T E S 


##installed py3 in venv, then works in sourced venv ***under py2.7***
## = = = = 

import pandas as pd
import pandas.io.api as pio

import skbio.stats.composition
from skbio.stats.composition import multiplicative_replacement as mult_rep
from skbio.stats.composition import ancom

#path = input('absolute path to data table (or press \'enter\' to dismiss):')
#if len(path) == 0:
#    data = '~/data/data100_otu.txt' 
#else:
#    data = path
#    # insert test to check that data is real / loadable

#data = '/home/eri/Downloads/otu100.txt'
#data = '/home/jfg/Dropbox/SilentGeno/R/R_RY/ry100_otu.txt'
data = '/home/eri/Downloads/test/ry100_otu.txt'

# input and test
data_otu = pio.read_table(data)
#data_otu = pd.DataFrame.transpose(data_otu)


# multiplicative replacement - adaptive pseudocounts for 'zero' observations
data_otu_rep = mult_rep(data_otu)
index = data_otu.index.values
columns = data_otu.columns.values
data_otu_rep = pd.DataFrame(data_otu_rep,index=index,columns=columns)

#data - 3 group: 55C/65C
groupA=pd.Series([
		'1','1','1','1','1','1','1','1','1','1',
		'1','1','1','1','1','1','1','1','1','1',
		'1','1','1','1','1','1','1','1','1','1',
		'1','1','1','1','1','1','1','1','1','1',
		'1','1','1','1','1','1','1','1','1','1',
		'1','1','1','1','1','1','1','1','1','1',
		'1','1','1','1','1','1','1','1','1','1',
		'1','1','1','1','1','1','1','1','1','1',
		'1','1','1','1','1','1','1','1','1','1',
		'1','1','1','1',
		
		'2','2','2','2','2','2','2','2','2','2',
		'2','2','2','2','2','2','2','2','2','2',
		'2','2','2','2','2','2','2','2','2','2',
		'2','2','2','2','2','2','2','2','2','2',
		'2','2','2','2','2','2','2','2','2','2',
		'2','2','2','2','2','2','2','2','2','2',
		
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3','3','3','3','3','3','3',
		'3','3','3','3'],index=[index])

##data - 3 group: B/C/D
#groupB=pd.Series(['1','1','1','2','2','2','3','3','3'],index=[index])



# Analysis of Composition Test
resultsA = ancom(data_otu_rep,groupA)
resultsA_df, resultsA_perc = resultsA
resultsA_df.to_csv('~/Bioinformatics/ancom/resultsdata_df.csv')
resultsA_df.to_csv(r'~/Bioinformatics/ancom/resultsdata_df.tsv', sep=',', mode='a')


	# superbasic out method - but WORKS
	text_file = open("escaperope.txt", "w")
	text_file.write("resultsA_df")
	text_file.close()

# breakdown of outcomes
resultsA_df_summ = resultsA_df.groupby(['Reject null hypothesis','W'])
resultsA_df_summ.size()
resultsA_df_summ.size.to_csv('resultsdata_df_summadata.csv')
	# superbasic strikes again
	text_file = open("escaperope_summary.txt", "w")
	text_file.write("resultsA_df")
	text_file.close()


## Analysis of Composition Test
#resultsB = ancom(data_otu_rep,groupB)
#resultsB_df, resultsB_perc = resultsB
#resultsB_df.to_csv('resultsB_df.csv')
## breakdown of outcomes
#resultsB_df_summ = resultsB_df.groupby(['Reject null hypothesis','W'])
#resultsB_df_summ.size()


##
##for i in comparisons:
##    result = ancom(data_otu_rep,i)
##    result_df, result_perc = result
##    return i = [result_df, result_perc]
