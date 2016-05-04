#!/bin/bash
#VarScan pipeline
#Rahul K. Das
#03/12/2016
#CAUTION!notice merged bam because sometimes there is no merged bam file
#<usage> bash varscan_pipeline.sh <samplename> <normalbam> <tumorbam>
#for normal bam and tumor bam, provide only the relative paths w.r.t the sample's Normal/Tumor directory;
#e.g. N-2/PTRIM.bam or Merged/PTRIM.bam etc

sample_name="$1"
normal_bam="$2"
tumor_bam="$3"

echo "[$(date):] Execution started for $sample_name on Normal/$normal_bam and Tumor/$tumor_bam"

projdir="/rawdata/projects/Radiogenomics"
fastafile="/results/referenceLibrary/tmap-f3/hg19/hg19.fasta"
varscan="/rawdata/software/varscan/VarScan.v2.4.1.jar"
bamreadcount="/rawdata/software/bam-readcount/bin/bam-readcount"
currdir=$(pwd)
workdir="$projdir/$sample_name/VarScan_Analysis/"
#echo $workdir

#make working directory
if [ -d $workdir ];
then
	echo "VarScan analysis directory already exists; not overwriting"
else
	mkdir -p $workdir
fi

cd $workdir

#convert primer trimmed merged normal and tumor bam files to pileup
normal_bam="$projdir/$sample_name/Normal/$normal_bam"
tumor_bam="$projdir/$sample_name/Tumor/$tumor_bam"

if [ -f "normal.pileup" ];
then	
	echo "normal pileup already exists, not rerunning samtools mpileup"
else
	echo "[$(date):] Converting normal bam to pileup.."
	samtools mpileup -q 1 -f $fastafile $normal_bam > normal.pileup
fi

if [ -f "tumor.pileup" ];
then
	echo "tumor pileup already exists, not rerunning samtools mpileup"
else
	echo "[$(date):] Converting tumor bam to pileup.."
	samtools mpileup -q 1 -f $fastafile $tumor_bam > tumor.pileup
fi
echo "[$(date):] Finished converting bam files to pileup.."


#run VarScan on normal and tumor pileup
varscan_command="java -jar $varscan somatic normal.pileup tumor.pileup output --min-var-freq 0.10 --strand-filter 1 --validation 1"
echo "[$(date):] Started varscan somatic run with command: $varscan_command"
$varscan_command
echo "[$(date):] Finished varscan somatic run.."

echo "[$(date):] Applying varscan somatic filters.."
#run varscan's processSomatic filter on snp/indel outputs to filter hc
java -jar $varscan processSomatic output.snp
java -jar $varscan processSomatic output.indel

#make the region files based on hc snps/indels; this will be used for next step: bam-readcount
awk '{if (NR>1) print $1"\t"$2"\t"$2}' output.snp.Somatic.hc > region.snp.Somatic.hc
awk '{if (NR>1) print $1"\t"int($2)-10"\t"int($2)+10}' output.indel.Somatic.hc > region.indel.Somatic.hc

#run bam-readcount on the hc snp/indel positions of tumor bam file 
bamreadcount_options="-q 1 -b 10"
echo "Running bam-readcount with following options: $bamreadcount_options"
$bamreadcount $bamreadcount_options -f $fastafile -l region.snp.Somatic.hc $tumor_bam > tumor_bam_readcount_snp.txt
$bamreadcount $bamreadcount_options -f $fastafile -l region.indel.Somatic.hc $tumor_bam > tumor_bam_readcount_indel.txt


#run varscan's fpfilter on hc.snp file
java -jar $varscan fpfilter output.snp.Somatic.hc tumor_bam_readcount_snp.txt --output-file output.snp.Somatic.hc.fpfilter.pass
java -jar $varscan fpfilter output.indel.Somatic.hc tumor_bam_readcount_indel.txt --output-file output.indel.Somatic.hc.fpfilter.pass


#extract 30x covered variants
#awk '{if (NR<2 || (int($5)+int($6)>=30 && int($9)+int($10)>=30)) print}' <output.snp.Somatic.hc.fpfilter.pass > output.snp.Somatic.hc.fpfilter.pass_30x

cd $currdir
echo "[$(date):] Execution finished.."


