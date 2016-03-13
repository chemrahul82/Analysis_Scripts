#intersects a varscan output file with the annovar-annotated file and makes a combined output file
#Rahul K. Das
#03.11.2016
#<usage>: python varscan_annovar.py <varscan file> <annovar file> <output file>

import sys
import time
def varscan_annovar(varfile,annfile,outfile):
	headerlist=['Chr','Position','Ref','Alt','Gene','AA-Change','Normal_Depth','Normal_Ref_Depth','Normal_Alt_Depth',\
			'Normal MAF','Tumor_Depth','Tumor_Ref_Depth','Tumor_Alt_Depth','Tumor MAF','Somatic-p-value','Cosmic68','LJB23_SIFT_Pred',\
			'LJB23_Polyphen2_HDIV_Pred','LJB23_Polyphen2_HVAR_Pred','LJB23_LRT_Pred','LJB23_MutationTaster_Pred',\
			'LJB23_MutationAssessor_Pred','LJB23_FATHMM_Pred','LJB23_RadialSVM_Pred','LJB23_LR_Pred',]
	header='\t'.join(headerlist)
	#print header
	
	#store varscan info
	with open(varfile, 'r') as f1:
		next(f1)
		fields1=[line for line in f1];
	
	#store annovar info
	with open(annfile, 'r') as f1:
		next(f1)
		fields2=[line for line in f1];
	
	#intersect and print
	with open(outfile, 'w') as fo:
		fo.write(header+'\n')
		for iline in fields1:
			line1 = iline.split('\t') #varscan
			string1=line1[0]+'_'+line1[1]
			for jline in fields2:
				line2 = jline.split('\t') #annovar
				string2=line2[0]+'_'+line2[1]
				if string1 == string2:
					newline=line1[0:4]
					newline.append(line2[6])#gene
					newline.append(line2[9].split(',')[0].split(':')[4][2:]) #AA-change
					newline.append(str(int(line1[4])+int(line1[5]))) #total normal depth
					newline=newline+line1[4:7]
					newline.append(str(int(line1[8])+int(line1[9]))) #total tumor depth
					newline=newline+line1[8:11]
					newline.append(line1[14])
					newline.append(line2[13])
					newline.append(line2[16])
					newline.append(line2[18])
					newline.append(line2[20])
					newline.append(line2[23])
					newline.append(line2[26])
					newline.append(line2[29])
					newline.append(line2[32])
					newline.append(line2[35])
					newline.append(line2[37])
					#print newline	
					fo.write('\t'.join(newline) + '\n')
					#print newline
					#time.sleep(30)

if __name__ == '__main__':
	f1=sys.argv[1]
	f2=sys.argv[2]
	f3=sys.argv[3]
	varscan_annovar(f1,f2,f3)
