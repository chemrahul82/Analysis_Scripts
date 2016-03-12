#!/bin/bash
#VarScan pipeline
#Rahul K. Das
#03/12/2016
#CAUTION!notice merged bam because sometimes there is no merged bam file
#<usage> bash varscan_pipeline.sh <samplename> <normalbam> <tumorbam>

sample_name=$1
normal_bam=$2
tumor_bam=$3
echo "[$(date):] Execution started for "$sample_name "on" $normal_bam "and" $tumor_bam

fastafile="/results/referenceLibrary/tmap-f3/hg19/hg19.fasta"
varscan="/rawdata/software/varscan/VarScan.v2.4.1.jar"
bamreadcount="/rawdata/software/bam-readcount/bin/bam-readcount"
currdir=$(pwd)
vardir="VarScan_Analysis"
mkdir -p "/rawdata/projects/Radiogenomics/$sample_name/$vardir"


echo "[$(date):] Converting normal bam to pileup.."
#convert primer trimmed merged normal and tumor bam files to pileup
samtools mpileup -q 1 -f $fastafile $normal_bam > normal.pileup
echo "[$(date):] Converting tumor bam to pileup.."
samtools mpileup -q 1 -f $fastafile $tumor_bam > tumor.pileup

echo "[$(date):] Started varscan somatic run.."
#run VarScan on normal and tumor pileup
java -jar $varscan somatic normal.pileup tumor.pileup output --strand-filter 1 --validation 1
echo "[$(date):] Finished varscan somatic run.."

echo "[$(date):] Applying varscan somatic filters.."
#run varscan's processSomatic filter on snp/indel outputs to filter hc
java -jar $varscan processSomatic output.snp.Somatic
java -jar $varscan processSomatic output.indel.Somatic

#make the region files based on hc snps/indels; this will be used for next step: bam-readcount
awk '{if (NR>1) print $1"\t"$2"\t"$2}' output.snp.Somatic.hc > region.snp.Somatic.hc
awk '{if (NR>1) print $1"\t"int($2)-10"\t"int($2)+10}' output.indel.Somatic.hc > region.indel.Somatic.hc

#run bam-readcount on the hc snp/indel positions of tumor bam file 
bamreadcount -q 1 -b 10 -f $fastafile -l region.snp.Somatic.hc $tumor_bam > tumor_bam_readcount_snp.txt
bamreadcount -q 1 -b 10 -f $fastafile -l region.indel.Somatic.hc $tumor_bam > tumor_bam_readcount_indel.txt


#run varscan's fpfilter on hc.snp file
java -jar $varscan fpfilter output.snp.Somatic.hc tumor_bam_readcount_snp.txt --output-file output.snp.Somatic.hc.fpfilter.pass
java -jar $varscan fpfilter output.indel.Somatic.hc tumor_bam_readcount_indel.txt --output-file output.indel.Somatic.hc.fpfilter.pass


#extract 30x covered variants
#awk '{if (NR<2 || (int($5)+int($6)>=30 && int($9)+int($10)>=30)) print}' <output.snp.Somatic.hc.fpfilter.pass > output.snp.Somatic.hc.fpfilter.pass_30x

echo "[$(date):] Execution finished.."
