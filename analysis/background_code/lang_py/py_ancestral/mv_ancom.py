##   N O T E S  -  
# see https://github.com/biocore/scikit-bio/blob/master/skbio/stats/composition.py
# for explanation and workthrough, especially use of the resultsX output. Briefly:
#
#	- resultsX_df holds hypothesis test, and W (number of SVs differed from).
#	  this shows occurrence of DA, and severity but not direction
#
#	- resultsX_perc holds percentile abundnances and therefore direction for 
#	  DA.


# pd.Series(...index=...) will be 'index' / 'column' depending when df transposed 
# surely improve on this: half the codework?

# idea to customise for each of the different projects: AJ, DW, MV. Feel free when free!
# thought could use kruskal as test for multiuple comparisons, but keeps erroring
#resultsB = ancom(sv_rep, groupB, significance_test = scipy.stats.kruskal)  # first failboat method 
#resultsB = ancom(sv_rep, groupB, significance_test = kruskal)		# second failboat method
# = = = = 

import pandas as pd
import pandas.io.api as pio
import skbio.stats.composition
from skbio.stats.composition import multiplicative_replacement as mult_rep
from skbio.stats.composition import ancom
# the piece we THOUGHT we were missing for multiple groups, but errors at 'all values equal'
from scipy.stats import kruskal
from scipy.stats import f_oneway  # not necess to specify at call, being explicit for clarity
# = = = = 

data = '/media/cart/Dropbox/SilentGeno/R/mv_r/mv_sv.txt'

# input and test
sv = pio.read_table(data)
#sv = pd.DataFrame.transpose(sv)    neutered for MV orientation

# multiplicative replacement - adaptive pseudocounts for 'zero' observations
sv_rep = mult_rep(sv)
sv_index = sv.index.values
sv_columns = sv.columns.values
sv_rep = pd.DataFrame(sv_rep,index=sv_index,columns=sv_columns)

# define groups for testing, as per sorting on the exported SV-TABLE
# MV - 2 group: ES/IS
groupA=pd.Series(['1','1','1','1','1','1','1','1','1','1','1','1','1','1','1',
		'2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2'],index=[sv_index])
# MV - 4 group: BES/CES/ISA/ISB
groupB=pd.Series(['1','1','1','1','1','1','2','2','2','2','2','2','2','2','2',
		'3','3','3','3','3','3','3','3','3','4','4','4','4','4','4','4','4','4',],index=[sv_index])


## Analysis of Composition Test - 2 Groups
resultsA = ancom(sv_rep,groupA)
resultsA_df, resultsA_perc = resultsA
resultsA_df.to_csv('resultsA_df.csv')  # writeout 
resultsA_perc.to_csv('resultsA_perc.csv')

# breakdown of outcomes, not an output in itself
resultsA_df_summ = resultsA_df.groupby(['W','Reject null hypothesis'])
resultsA_df_summ.size()

# once ID'd as DA, investigate direction of DA at percentile 50 (i.e. median value) like so:
resultsA_perc[50.0].loc['mvsv0001']
# or in more detail by indexing the SV directly for all abundance percentiles in all groups:
resultsA_perc.loc['mvsv0001']



## Analysis of Composition Test - 4 Groups (i.e. >2 groups)
resultsB = ancom(sv_rep,groupB,significance_test=f_oneway)  		# new (explicit) failboat recommended by documentation
resultsB_df, resultsB_perc = resultsB
resultsB_df.to_csv('resultsB_df.csv')
resultsB_perc.to_csv('resultsB_perc.csv')

resultsB_df_summ = resultsB_df.groupby(['W','Reject null hypothesis'])
resultsB_df_summ.size()

resultsB_perc[50.0].loc['mvsv0001']
resultsB_perc.loc['mvsv0001']


## handy graphics command? boxplots in scikit bio?
