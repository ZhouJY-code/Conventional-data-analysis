
rawdata="6h-wt-1_FRAS210206459-2r"
data="wt_6h_1"
ref="single-transcript-info"
ref_G="Danio_rerio.Zv9.79.gtf"

export PERL5LIB=/asnas/yangyg_group/wangshzh/software/czplib-master:/software/biosoft/software/bioperl:/asnas/yangyg_group/wangshzh/software/czplib-master/lib/perl5/x86_64-linux-thread-multi/:/software/biosoft/software/bioperl/lib/perl5
export PATH=$PATH:/software/biosoft/software/python/python2.7_2018_12/bin/:/software/biosoft/software/bwa-0.7.17:/software/biosoft/software/bedtools2.28/bedtools2/bin:/software/biosoft/software/fastqc/FastQC:/software/biosoft/software/tophat-2.1.1.Linux_x86_64:/software/biosoft/software/hisat2-2.0.5


cd ${data}


####quality control
mkdir clean_data
cd clean_data
mkdir fastqc

fastqc $rawdata/*_1.fq.gz -t 4 -o ./fastqc
fastqc $rawdata/*_2.fq.gz -t 4 -o ./fastqc

###remove adapter
cutadapt -a GATCGGAAGA -A GATCGGAAGA -o cutadapt_1.fastq.gz -p cutadapt_2.fastq.gz $rawdata/*_1.fq.gz $rawdata/*_2.fq.gz

fastqc cutadapt_1.fastq.gz -t 4 -o ./fastqc
fastqc cutadapt_2.fastq.gz -t 4 -o ./fastqc

###trim low quality reads，quality 20， min length 18 
java -Xmx4g -jar /software/biosoft/software/Trimmomatic-0.36/trimmomatic-0.36.jar  PE -phred33 ./cutadapt_1.fastq.gz ./cutadapt_2.fastq.gz ./trim1_1.fastq.gz ./unpaired1_1.fastq.gz ./trim1_2.fastq.gz ./unpaired1_2.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:18


fastqc trim1_1.fastq.gz -t 4 -o ./fastqc
fastqc trim1_2.fastq.gz -t 4 -o ./fastqc
rm cutadapt_1.fastq.gz cutadapt_2.fastq.gz
rm unpaired1_1.fastq.gz unpaired1_2.fastq.gz



###########Genome mapping
cd ${data}
mkdir hisat2-GRCh38
hisat2 -p 4 --dta  -x /pnas/yangyg_group/yangxin/reference/zebrafish/genome.fa -1 ./clean_data/trim1_1.fastq.gz -2 ./clean_data/trim1_2.fastq.gz -S ./hisat2-GRCh38/hisat2.sam

cd hisat2-GRCh38

###Read counts summary
echo "reads with no -q:" > readsCount
grep -v "^@" hisat2.sam|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
###q20
samtools view -hb hisat2.sam > hisat2.bam
samtools view hisat2.bam -q 20 -h > hisat2_20.sam


###unique mapping reads
grep -E "^@|NH:i:1$|NH:i:1[^0-9]" hisat2_20.sam|awk '$1~/^@/ || $1~/^A0/' > uniqmap.sam
echo "reads with -q:" >> readsCount
grep -v "^@" hisat2_20.sam|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
rm hisat2_20.sam
rm hisat2.sam

echo "reads with uniqmap:" >> readsCount
cut -f1 uniqmap.sam|awk '!x[$0]++'|wc -l >> readsCount
echo "pair reads with uniqmap:" >> readsCount
awk '$7~/=/' uniqmap.sam|cut -f1|awk '!x[$0]++'|wc -l >> readsCount



samtools view -S uniqmap.sam -b -@ 4 -o uniqmap.bam
bedtools bamtobed -i uniqmap.bam -split|awk '{print"chr"$0}' > uniqmap.bed

rm uniqmap.sam


intersectBed -a uniqmap.bed -b /pnas/yangyg_group/yangxin/reference/zebrafish/single-transcript-info -wa -wb -f 0.51|awk '$6==$10' > tmp.SG2


echo "2h-dmso-1 ncRNA+mRNA reads:" >> readsCount
awk '!a[$4]++' tmp.SG2 |wc >> readsCount
echo "2h-dmso-1 mRNA reads:" >> readsCount
awk '$13=="protein_coding" && !a[$4]++' tmp.SG2|wc >> readsCount
echo "2h-dmso-1 ncRNA reads:" >> readsCount
awk '$13!="protein_coding" && !a[$4]++' tmp.SG2|wc >> readsCount

echo "2h-dmso-1 mRNA region reads:" >> readsCount
cat tmp.SG2|awk '$13=="protein_coding"'|cut -f14|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount
rm tmp.SG2
rm uniqmap.bed

##htSeq
samtools sort -n uniqmap.bam -@ 4 -o sort.bam

mkdir ./htSeq
htseq-count -m union -f bam -s no sort.bam ${ref_G} > ./htSeq/union_no.out
rm sort.*



#######featureCounts
samtools sort uniqmap.bam -@ 4 -o sort.bam
samtools index sort.bam
rm uniqmap.bam

mkdir ./featurecount
/xtdisk/yangyg_ncov/gaochchA/SOFTWARE/subread-1.6.5-source/bin/featureCounts -T 4 -p -a $ref_G  -t exon -g gene_id -o ./featurecount/counts.txt sort.bam

