#compare varscan predicted final TP nonsyn/stopgain snps between QC-merged and all-merged runs
#usage: python qc_all_compare.py <annovar annotated file for QC-merged> <annovar annotated file for all-merged> <merged fpfilterpass> <allmerged fpfilterpass>
#variants that are nonsyn/stopgain and no gemline evidence in databases

import sys
def qc_all_compare(qcannfile,nonqcannfile,qcvarfile,nonqcvarfile):
	with open(qcannfile, 'r') as f1:
		next(f1)
		string1=[line.split('\t')[0]+'_'+line.split('\t')[1] for line in f1]
	with open(nonqcannfile, 'r') as f2:  
		next(f2)
		string2=[line.split('\t')[0]+'_'+line.split('\t')[1] for line in f2]

	#intersection
	common=(list(set(string1) & set(string2)))
	#unique to QC-merged
	unique2qc=(list(set(string1)-set(string2)))
	#unique to all-merged
	unique2all=(list(set(string2)-set(string1)))
	
	#print info for Venn Diagram
	print len(unique2qc), len(common), len(unique2all) 

	fo1=open('unique2qc.txt', 'w')
	fo2=open('unique2all.txt', 'w')
	fo3=open('common.txt', 'w')
	
	#write varscan output for variants that are unique to QC-merged
	with open(qcvarfile, 'r') as f3:
		lines=[line for line in f3]
	
	for idx, iline in enumerate(lines): 
	 	if idx == 0:
			fo1.write(iline)	
		for entries in unique2qc:		
			if entries==iline.split('\t')[0]+'_'+iline.split('\t')[1]:
				fo1.write(iline)	
	
	#write varscan output for variants that are unique to all-merged
	with open(nonqcvarfile, 'r') as f4:
		lines=[line for line in f4]
	
	for idx, iline in enumerate(lines): 
	 	if idx == 0:
			fo2.write(iline)	
		for entries in unique2all:		
			if entries==iline.split('\t')[0]+'_'+iline.split('\t')[1]:
				fo2.write(iline)	
	

	#write varscan output for variants that are common
	with open(qcvarfile, 'r') as f4:
		lines=[line for line in f4]
	
	for idx, iline in enumerate(lines): 
	 	if idx == 0:
			fo3.write(iline)	
		for entries in common:		
			if entries==iline.split('\t')[0]+'_'+iline.split('\t')[1]:
				fo3.write(iline)	
	


	#return common


if __name__ == '__main__':
	f1=sys.argv[1]
	f2=sys.argv[2]
	f3=sys.argv[3]
	f4=sys.argv[4]
	qc_all_compare(f1,f2,f3,f4)
