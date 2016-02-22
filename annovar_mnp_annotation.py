#!/usr/bin/env python2.7

"""
function for annotating AA-changes for MNP variants; runs Annovar a second time using coding_change.pl script
Author: Rahul K. Das
Date: Feb 17, 2016

"""

import os, sys
def mnp_aachange(sv_dir,software_dir):

	#prefix for annovar output files
	#prfx='temppp'
	prfx = 'annovar_somatic_table'
	#exonic variant function file created by 1st run of annovar
	exnvarfxnfile = '%s/%s.refGene.exonic_variant_function' %(sv_dir, prfx)

	#final annovar multianno.vcf file that contains annotations for mnps after the 2nd run
	multiannovcffinal = '%s/%s.hg19_multianno_final.vcf' %(sv_dir, prfx)	

	#two temporary files 
	tmpvarfxnfile = '%s/%s.temp.refGene.exonic_variant_function' %(sv_dir, prfx)
	tmpmnpannotfile = '%s/mnp.log' %(sv_dir)

	#command for running coding_change.pl
	mnpannotcomm = '%s/annovar/coding_change.pl \
	%s \
	/rawdata/software/annovar/humandb_ucsc/hg19_refGene.txt \
	/rawdata/software/annovar/humandb_ucsc/hg19_refGeneMrna.fa > %s' %(software_dir, tmpvarfxnfile, tmpmnpannotfile)

	mf = open(multiannovcffinal, 'w') #initialize the file handle for final vcf output
	#scan through multianno output vcf from 1st run of annovar
	with open('%s/%s.hg19_multianno.vcf' %(sv_dir, prfx), 'r') as annovcf:
		for idx,line in enumerate(annovcf):
			fields=(line.split('\t'))
			#annotations of 8th field in line
			annotations=(fields[7].split(';'))
			found_mnp = False

			#we only only run coding_change.pl for exonic mnp variants
			if (("Func.refGene=exonic" in annotations) and ("TYPE=mnp" in annotations) \
				and ("AAChange.refGene=UNKNOWN" not in annotations)):
				#search for TYPE=MNP"; we need to rerun annovar on these
				found_mnp=True
				#line number for the mnp variant	
				lineidx = 'line'+ '%i' %(idx+1)
					
				#now grab the line for the mnp variant from annovar exonic variant function file of 1st run and create a temp file		
				with open(exnvarfxnfile, 'r') as fi:	
					for line in fi:
						if line.split('\t')[0]== lineidx:
							with open(tmpvarfxnfile, 'w') as fo:
								fo.write(line)

				#now run coding_change.pl to get the AA-change annotations for MNP and save the output to a temp file	
				os.system(mnpannotcomm)
				#extract the transcript and AA-change info from the temp file
				with open(tmpmnpannotfile, 'r') as fi:
					for line in fi:
						if "position" in line:
							attributes=line.split(' ')
							mnp_exon_fxn = attributes[3] #detailed exonic function
							mnp_transcript = attributes[1]
							#the AA change; Annovar gives an annotation that is a sentence long; trim it	
							mnp_AAChange = [i+j+k for i,j,k in zip(attributes[9],attributes[6].split('-')[0:],attributes[11])]
							#mnp_AAChange = ' '.join(attributes[5:]).strip('\n')
				#remove the two temp files	
				os.remove(tmpvarfxnfile)	
				os.remove(tmpmnpannotfile)	
			
			if found_mnp==True:
				for ia, annot in enumerate(annotations):
					if "ExonicFunc" in annot:
						#append the detailed exonic function annotation to the existing one
						annotations[ia]=annot+'_'+mnp_exon_fxn
					if "AAChange." in annot:
						#change the annotation for the AAChange.refGene field of the mnp
						newstring=annot.split(',')[0].split(':')[0]+':'+mnp_transcript+':'+\
						annot.split(',')[0].split(':')[2]+':'+annot.split(',')[0].split(':')[3]+':'+'p.'+' & '.join(map(str,mnp_AAChange))
						#update the list with new annotation
						annotations[ia] = newstring
												
							
			#update the 8th field of the line for mnp
			fields[7]=';'.join(annotations)
			#with open(multiannovcffinal, 'w') as mf:
			mf.write('\t'.join(fields))				
	mf.close()
	return;
if __name__ == "__main__":
	#how to run
	#annovar_mnp_annotation.py <somatic_var_dir> <annovar_dir>
	somaticvar_dir = sys.argv[1]
	annovar_dir = sys.argv[2]
	mnp_aachange(somaticvar_dir, annovar_dir)					
		
