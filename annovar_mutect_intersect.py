#intersects a mutect output file with the annovar-annotated file and makes a combined output file
#Rahul K. Das
#03.11.2016
#<usage>: python mutect_annovar.py <mutect file> <annovar file> <output file>

import sys
import time
def mutect_annovar(mutfile,annfile,outfile):
	headerlist=['Chr','Position','Context','Ref','Alt','Gene','AA-Change','Normal_Depth','Normal_q20','Normal_Ref_Depth','Normal_Alt_Depth',\
			'Normal_MAF','Tumor_Depth','Tumor_q20','Tumor_Ref_Depth','Tumor_Alt_Depth','Tumor_MAF','Power','Normal_Power','Tumor_Power',\
			'Cosmic68','LJB23_SIFT_Pred','LJB23_Polyphen2_HDIV_Pred','LJB23_Polyphen2_HVAR_Pred','LJB23_LRT_Pred','LJB23_MutationTaster_Pred',\
			'LJB23_MutationAssessor_Pred','LJB23_FATHMM_Pred','LJB23_RadialSVM_Pred','LJB23_LR_Pred',]
	header='\t'.join(headerlist)
	#print header
	
	#store mutect info
	with open(mutfile, 'r') as f1:
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
			line1 = iline.split('\t') #mutect
			string1=line1[0]+'_'+line1[1]
			for jline in fields2:
				line2 = jline.split('\t') #annovar
				string2=line2[0]+'_'+line2[1]
				if string1 == string2:
					newline=line1[0:5]
					newline.append(line2[6])#gene
					newline.append(line2[9].split(',')[0].split(':')[4][2:]) #AA-change
					newline.append(str(int(line1[37])+int(line1[38]))) #total normal depth
					newline.append(line1[36]) #n_q20
					newline=newline+line1[37:39]
					newline.append(line1[35]) #n_maf
					newline.append(str(int(line1[25])+int(line1[26]))) #total tumor depth
					newline.append(line1[24]) #t_q20
					newline=newline+line1[25:27]
					newline.append(line1[21]) #t_maf
					newline.append(line1[10]) #power
					newline.append(line1[12]) #n_power
					newline.append(line1[11]) #t_power
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
	mutect_annovar(f1,f2,f3)
