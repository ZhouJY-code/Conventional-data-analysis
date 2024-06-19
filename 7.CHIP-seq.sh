
##preprocess
path="CHIP"
i="Rep1"

cd $path

trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o ./clean $i_1.fastq.gz $i_2.fastq.gz

cd $path/$i/clean
fastqc -o $Fpath/$i/qc/ $i_1_val_1.fq.gz $i_2_val_2.fq.gz

##Mapping
cd $path/clean

bwa mem -t 30 /home/Reference/T2T_V1.1-CoV/bwa.index/T2T-CoV2.fa $i_1_val_1.fq.gz $i_2_val_2.fq.gz > $path/$i/align/$i.sam

cd $path/$i/align

samtools sort -@ 15 -m 1G $i.sam -o $i.sortP.bam
samtools index $i.sortP.bam

##Filter
cd $path/$i/align

sambamba markdup  --overflow-list-size 600000 -r $i.sortP.bam  -t 16 $i.rmdup.bam
samtools view  -h  -f 2 -q 30 $i.rmdup.bam |grep -v chrM |samtools sort  -O bam  -@ 15 -o - >$i.last.bam
samtools index $i.last.bam
samtools flagstat $i.last.bam > bam_last.stat
         
##CallPeak
macs2 callpeak -t $path/$i/align/$i.last.bam -c $path/$i/align/bwa.T2T.input.last.bam -f BAM -g 3.055e9 -n H3K27ac -B -q 0.05
