
path="ATAC-seq"
i="Rep1"

#Quality control
cd $path/Cleandata && mkdir fastQC

fastp -i $i_R1.fq.gz -o $i_R1_fastp.fq.gz -I $i_R2.fq.gz  -O  $i_R2_fastp.fq.gz  --thread=16  --compression=4
fastqc -t 16 -o ./fastQC $i_R1_fastp.fq.gz $i_R2_fastp.fq.gz 		

#Mapping
mkdir bwa-T2T && cd bwa-T2T

bwa mem  -t 16 ./bwa.index/T2T-CoV2.fa $path/Cleandata/$i/$i_R1_fastp.fq.gz $path/Cleandata/$i/$i_R2_fastp.fq.gz  > $i.sam

#File process
samtools sort -@ 16 -m 1G $i.sam -o $i.sortP.bam
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true  I=$i.sortP.bam O=$i.sortP_rmdupPICARD.bam M=$i.marked_dup_metrics.txt
	
samtools index $i.sortP_rmdupPICARD.bam
samtools view -@ 10 -h $i.sortP_rmdupPICARD.bam |grep -v -e 'XA:Z:' -e 'SA:Z:' |samtools view -@ 10 -b > $i.sortP_rmdupPICARD_uniq.bam	

#Call peak
mkdir MACS2 && cd MACS2
samtools sort -@ 16 -m 1G ../bwa-T2T/$i.sam -o $i.sortP.bam
macs2 callpeak -t $path/picard/$i.sortP_rmdupPICARD.bam -n $i --shift -100 --extsize 200 --nomodel -B --SPMR -g 3.055e9	


