#performs a two-sample t-test to determine the statistical significance of the difference in genotypes at a given position 
#between two runs of same sample
# if the p-value < alpha, the two difference in the genotype at that position is not counted as error
#Rahul K. Das
#<Usage>: python cal_3x3_error_p_val_adjusted.py <project directory> <sampleName> < Run 1> <Run 2> <RunType1> <RunType2> <alpha> 

from __future__ import division
import os, sys
from scipy import stats
import json

def cal_pval_3x3_error(projDir, sampleName, run1, run2, runType1, runType2, alpha):
	
	#counters for mismatch
	WT_HET = 0
	WT_HOM = 0
	HET_WT = 0
	HET_HOM = 0
	HOM_WT = 0
	HOM_HET = 0

	matchVarFile=os.path.join(projDir, sampleName, 'QC', 'chr1'+run1+'vs'+run2,'matched_variantschr1.csv')
	jsonFile=os.path.join(projDir, sampleName, 'QC','results_QC.json')
	
        with open(matchVarFile, 'r') as f:
		next(f)
		for line in f:
			fields = line.split('\t')
			#exclude the positions with no GT assigments
			if fields[4] != '.' and fields[8] != '.':
				minor1 = int(fields[6]) #alt reads for run1
				ref1 = int(fields[7]) #ref reads for run1
				minor2 = int(fields[10]) #alt reads for run2
				ref2 = int(fields[11]) #ref reads for run2
				list1 = ([1]*minor1) + ([0]*ref1) #allele distribution for run2
				list2 = ([1]*minor2) + ([0]*ref2) #allele distribution for run2
				
                                if fields[4] == 'WT' and fields[8] == 'HET':
                                    p_val = stats.ttest_ind(list1,list2)[1] # t-test; assuming equal variance
                                    if float(p_val) < float(alpha):
                                        WT_HET += 1
                                if fields[4] == 'WT' and fields[8] == 'HOM':
                                    p_val = stats.ttest_ind(list1,list2)[1] # t-test; assuming equal variance
                                    if float(p_val) < float(alpha):
                                        WT_HOM += 1
                                if fields[4] == 'HET' and fields[8] == 'WT':
                                    p_val = stats.ttest_ind(list1,list2)[1] # t-test; assuming equal variance
                                    if float(p_val) < float(alpha):
                                        HET_WT += 1
                                if fields[4] == 'HET' and fields[8] == 'HOM':
                                    p_val = stats.ttest_ind(list1,list2)[1] # t-test; assuming equal variance
                                    if float(p_val) < float(alpha):
                                        HET_HOM += 1
                                if fields[4] == 'HOM' and fields[8] == 'WT':
                                    p_val = stats.ttest_ind(list1,list2)[1] # t-test; assuming equal variance
                                    if float(p_val) < float(alpha):
                                        HOM_WT += 1
                                if fields[4] == 'HOM' and fields[8] == 'HET':
                                    p_val = stats.ttest_ind(list1,list2)[1] # t-test; assuming equal variance
                                    if float(p_val) < float(alpha):
                                        HOM_HET += 1
	
        
        with open(jsonFile) as data_file:    
                data = json.load(data_file)
        total_bases = int(data["QC_comparisons"]["chr1"]["%s_%s" %(runType1,runType2)]["%svs%s" %(run1,run2)]["total_eligible_bases"])        
        
        if runType1 == runType2:
            total_errors = WT_HET+WT_HOM+HET_WT+HET_HOM+HOM_WT+HOM_HET
        else:
            total_errors = HET_WT+HOM_WT+HOM_HET
        
        error_rate = total_errors/total_bases

	print WT_HET, WT_HOM, HET_WT, HET_HOM, HOM_WT, HOM_HET
        #print total_bases, total_errors, error_rate
        print 'total errors is %d: ' %total_errors +'error rate is %s: ' %error_rate +'at significance level %s' %alpha
				
if __name__ == "__main__":
	projDir = sys.argv[1]
	sampleName = sys.argv[2]
	run1 = sys.argv[3]
	run2 = sys.argv[4]
        runType1 = sys.argv[5]
        runType2 = sys.argv[6]
	alpha = sys.argv[7]
	cal_pval_3x3_error(projDir, sampleName, run1, run2, runType1, runType2, alpha) 


