

export PERL5LIB=/asnas/yangyg_group/wangshzh/software/czplib-master:/software/biosoft/software/bioperl:/asnas/yangyg_group/wangshzh/software/czplib-master/lib/perl5/x86_64-linux-thread-multi/:/software/biosoft/software/bioperl/lib/perl5
export PATH=$PATH:/software/biosoft/software/python/python2.7_2018_12/bin/:/software/biosoft/software/bwa-0.7.17:/software/biosoft/software/bedtools2.28/bedtools2/bin:/software/biosoft/software/fastqc/FastQC:/software/biosoft/software/tophat-2.1.1.Linux_x86_64:/software/biosoft/software/hisat2-2.0.5

data=iCLIP

R1=HJ-4_FKDL210011848-1a_1.fq.gz
thread=6


cd ${data}
mkdir fastqc

fastqc $R1 -t $thread -o ./fastqc

#####Trimming of 3' linker sequences
cutadapt -a AGATCGGAAGAG --discard-untrimmed -o S1_cutadapt_1.fq.gz $R1

fastqc S1_cutadapt_1.fq.gz -t $thread -o ./fastqc

cutadapt -a "A{10}" -u 3 -e 0.1 -m 18 -M 60 S1_cutadapt_1.fq.gz  -o S2_cutadapt_1.fq.gz

fastqc S2_cutadapt_1.fq.gz -t $thread -o ./fastqc

#Quality control of fastq file
java -Xmx4g -jar /software/biosoft/software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 -threads $thread S2_cutadapt_1.fq.gz S3_zebrafish-trim_1.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:18
fastqc S3_zebrafish-trim_1.fq.gz -t $thread -o ./fastqc

GTF=/pnas/yangyg_group/yangxin/reference/zebrafish/Danio_rerio.Zv9.79.gtf
genome=/pnas/yangyg_group/yangxin/reference/zebrafish/BowtieIndex/genome.fa
ref=/pnas/yangyg_group/yangxin/reference/zebrafish/single-transcript-info-mRNA-nochr



######single reads mapping
mkdir tophat
cd tophat
tophat --bowtie1 -p $thread -G $GTF -o ./ $genome ../S3_zebrafish-trim_1.fq.gz

samtools view accepted_hits.bam -q 20 -h > uniqmap.sam     #-q 20
samtools view -S uniqmap.sam -b -o uniqmap.bam
bedtools bamtobed -i uniqmap.bam -split > uniqmap.bed
#rm accepted_hits.sam
#rm uniqmap.sam


samtools sort -n uniqmap.bam -o accepted_hits_unique_name_sorted.bam

#exon reads
htseq-count -f bam accepted_hits_unique_name_sorted.bam -s no -m union -t exon -r name -a 20 ${GTF} > zebarfish-Single-htseqCount.txt

echo "reads with no -q:" > readsCount-genecode
samtools view accepted_hits.bam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount-genecode

samtools view accepted_hits.bam -q 20 -h > hisat2_20.sam
echo "reads with -q:" >> readsCount-genecode
cat hisat2_20.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount-genecode
rm hisat2_20.sam

echo "reads with uniqmap:" >> readsCount-genecode
samtools view accepted_hits_unique_name_sorted.bam|cut -f1|awk '!x[$0]++'|wc -l >> readsCount-genecode
echo "pair reads with uniqmap:" >> readsCount-genecode
samtools view accepted_hits_unique_name_sorted.bam|awk '$7~/=/'|cut -f1|awk '!x[$0]++'|wc -l >> readsCount-genecode


bedtools bamtobed -i accepted_hits_unique_name_sorted.bam -split > uniqmap.bed
intersectBed -a uniqmap.bed -b ${ref} -wa -wb -f 0.51 > tmp.SG2
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode

echo "zebarfish mapping reads:" >> readsCount-genecode
cat uniqmap.bed|awk '!x[$4]++'|wc >> readsCount-genecode
echo "zebarfish ncRNA reads:" >> readsCount-genecode
cat tmp.SG2 |awk '$13 != "protein_coding" && !a[$4]++'|wc >> readsCount-genecode
echo "zebarfish mRNA reads:" >> readsCount-genecode
cat tmp.SG2 |awk '$13 == "protein_coding" && !a[$4]++'|wc >> readsCount-genecode


echo "zebarfish all region reads:" >> readsCount-genecode
cat tmp.SG2|cut -f14|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode

echo "zebarfish region reads only for protein_coding:" >> readsCount-genecode
cat tmp.SG2|awk '$13 == "protein_coding"'|cut -f14|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode

rm tmp.SG2
rm uniqmap.bed




