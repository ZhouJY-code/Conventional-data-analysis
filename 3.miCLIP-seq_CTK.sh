
export LD_LIBRARY_PATH=/software/biosoft/software/zlib1.2.7/lib:/software/biosoft/software/gcc4.9_install/gcc4.9.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/software/biosoft/software/tbb44_20160316oss/build/linux_intel64_gcc_cc4.4.6_libc2.12_kernel2.6.32_release:/software/biosoft/software/zlib1.2.7/lib:/software/biosoft/software/gcc4.9_install/gcc4.9.2/lib64:$LD_LIBRARY_PATH
export PERL5LIB=/software/biosoft/bioperl/share/perl5:/share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/Perl_Lib/lib64/perl5:/share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/Perl_Lib/CIMS/czplib:$PERL5LIB

##Preprocessing process
##Keep read length within a short fragment length
cd /share_bio/unisvx1/yangyg_group/chenysh/ZhouQ_TESIS-m6A/NEW_ANALYSIS/P12-WT-miCLIP-Total

mkdir script_stout
mkdir clean_data
mkdir ./clean_data/fastqc
cd clean_data


gunzip -c /share_bio/disk8/yangyg_group/sunbf/ZhouQ_TESTIS-m6A/Sample_P12-WT-miCLIP/*_R1_001.fastq.gz > P12-WT-miCLIP_1.fastq
gunzip -c /share_bio/disk8/yangyg_group/sunbf/ZhouQ_TESTIS-m6A/Sample_P12-WT-miCLIP/*_R2_001.fastq.gz > P12-WT-miCLIP_2.fastq

fastqc P12-WT-miCLIP_1.fastq -t 4 -o ./fastqc
fastqc P12-WT-miCLIP_2.fastq -t 4 -o ./fastqc


#########################################################
##Trimming of 3' linker sequences
##Remove adaptors, and quality control (standards are very low)
##barcode 9nt+tag 15nt

fastx_clipper -a AGATCGGAAGAGCACACG -l 24 -n -i P12-WT-miCLIP_1.fastq -Q 33|fastq_quality_trimmer -t 5 -l 24 -Q 33 -o trim_P12-WT-miCLIP_1.fastq
fastx_clipper -a AGATCGGAAGAGCGTCGT -l 24 -n -i P12-WT-miCLIP_2.fastq -Q 33|fastq_quality_trimmer -t 5 -l 24 -Q 33 -o trim_P12-WT-miCLIP_2.fastq
rm P12-WT-miCLIP_1.fastq P12-WT-miCLIP_2.fastq

fastqc trim_P12-WT-miCLIP_1.fastq -t 4 -o ./fastqc
fastqc trim_P12-WT-miCLIP_2.fastq -t 4 -o ./fastqc


#########################################################
##Filter out the long fragments and keep the short ones
##The BARCODE conservative site at the end is required to be correct

awk 'BEGIN{permision=0}{if(NR%4==1){tmp=$0}else{if(NR%4==2){if(length($0)<=60 && $0~/^....TG/){print tmp; print $0; permission=1}else{permission=0}}else{if(permission==1){print $0}} }}' trim_P12-WT-miCLIP_1.fastq > filter_trim_P12-WT-miCLIP_1.fastq
awk 'BEGIN{permision=0}{if(NR%4==1){tmp=$0}else{if(NR%4==2){if(length($0)<=60 && $0~/CA....$/){print tmp; print $0; permission=1}else{permission=0}}else{if(permission==1){print $0}} }}' trim_P12-WT-miCLIP_2.fastq > filter_trim_P12-WT-miCLIP_2.fastq
gzip trim_P12-WT-miCLIP_1.fastq
gzip trim_P12-WT-miCLIP_2.fastq

fastqc filter_trim_P12-WT-miCLIP_1.fastq -t 4 -o ./fastqc
fastqc filter_trim_P12-WT-miCLIP_2.fastq -t 4 -o ./fastqc


#########################################################
##Collapse exact duplicates
##remove pcr duplication

perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/fastq2collapse.pl filter_trim_P12-WT-miCLIP_1.fastq deduplicate_P12-WT-miCLIP_1.fastq 
perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/fastq2collapse.pl filter_trim_P12-WT-miCLIP_2.fastq deduplicate_P12-WT-miCLIP_2.fastq 
gzip filter_trim_P12-WT-miCLIP_1.fastq
gzip filter_trim_P12-WT-miCLIP_2.fastq

fastqc deduplicate_P12-WT-miCLIP_1.fastq -t 4 -o ./fastqc
fastqc deduplicate_P12-WT-miCLIP_2.fastq -t 4 -o ./fastqc


#########################################################
##Strip barcode
##remove barcode

perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/stripBarcode.pl -format fastq -len 9 deduplicate_P12-WT-miCLIP_1.fastq temp_1.fastq
awk '{if(NR%4==1){split($1,a,"#");print a[1]"#1#"a[3]}else{print $0}}' temp_1.fastq > debarcode_P12-WT-miCLIP_1.fastq 
gzip deduplicate_P12-WT-miCLIP_1.fastq
rm temp_1.fastq
fastqc debarcode_P12-WT-miCLIP_1.fastq -t 4 -o ./fastqc

fastx_reverse_complement -Q 33 -i deduplicate_P12-WT-miCLIP_2.fastq -o deduplicate_P12-WT-miCLIP_2r.fastq
perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/stripBarcode.pl -format fastq -len 9 deduplicate_P12-WT-miCLIP_2r.fastq temp_2.fastq
awk '{if(NR%4==1){split($1,a,"#");print a[1]"#2#"a[3]}else{print $0}}' temp_2.fastq > debarcode_P12-WT-miCLIP_2r.fastq 
gzip deduplicate_P12-WT-miCLIP_2.fastq
gzip deduplicate_P12-WT-miCLIP_2r.fastq
rm temp_2.fastq
fastqc debarcode_P12-WT-miCLIP_2r.fastq -t 4 -o ./fastqc

cat debarcode_P12-WT-miCLIP_1.fastq debarcode_P12-WT-miCLIP_2r.fastq > debarcode_P12-WT-miCLIP_final.fastq
gzip debarcode_P12-WT-miCLIP_1.fastq
gzip debarcode_P12-WT-miCLIP_2r.fastq


#########################################################
##Read quality filtering

perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/fastq_filter.pl -v -if sanger -f mean:0-33:20 -of fastq debarcode_P12-WT-miCLIP_final.fastq qf_P12-WT-miCLIP_final.fastq 
rm debarcode_P12-WT-miCLIP_final.fastq
fastqc qf_P12-WT-miCLIP_final.fastq -t 4 -o ./fastqc




#########################################################
##Mapping

cd /share_bio/unisvx1/yangyg_group/chenysh/ZhouQ_TESIS-m6A/NEW_ANALYSIS/P12-WT-miCLIP-Total

##mapping to bwa
mkdir bwa
cd bwa

bwa aln -t 4 -n 0.06 -q 20 /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Mouse/ENSEMBLE_68/index/genome.fa ../clean_data/qf_P12-WT-miCLIP_final.fastq > P12-WT-miCLIP.sai
bwa samse /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Mouse/ENSEMBLE_68/index/genome.fa P12-WT-miCLIP.sai ../clean_data/qf_P12-WT-miCLIP_final.fastq > P12-WT-miCLIP.sam



samtools view P12-WT-miCLIP.sam -q 20 -h -S > P12-WT-miCLIP_q20.sam
grep -E "^@|XT:A:U" P12-WT-miCLIP_q20.sam|awk '$1~/^@/ || $1~/^K00386/' > uniqmap.sam

echo "reads with no -q:" > readsCount
cat P12-WT-miCLIP.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
samtools view -S P12-WT-miCLIP.sam -b -o P12-WT-miCLIP.bam

echo "reads with -q:" >> readsCount
cat P12-WT-miCLIP_q20.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
samtools view -S P12-WT-miCLIP_q20.sam -b -o P12-WT-miCLIP_q20.bam
rm P12-WT-miCLIP_q20.sam

echo "reads with uniqmap:" >> readsCount
cat uniqmap.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount

samtools view -S uniqmap.sam -b -o uniqmap.bam
bedtools bamtobed -i uniqmap.bam -split > uniqmap.bed


echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount


cat uniqmap.bed|awk -v OFS="\t" '{print $0}'|intersectBed -a - -b /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Mouse/ENSEMBLE_68/GeneList/single-transcript-chr -wa -wb -f 0.5 > tmp.P12-WT-miCLIP
echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount

echo "P12-WT-miCLIP mapping reads:" >> readsCount
cat uniqmap.bed|awk '!x[$4]++'|wc >> readsCount
echo "P12-WT-miCLIP ncRNA+mRNA reads:" >> readsCount
cat tmp.P12-WT-miCLIP |awk '!a[$4]++'|wc >> readsCount
echo "P12-WT-miCLIP mRNA region reads:" >> readsCount

cat tmp.P12-WT-miCLIP|awk '$13=="mRNA"'|cut -f14|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount
echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount
rm tmp.P12-WT-miCLIP



#########################################################
export PERL5LIB=/software/biosoft/bioperl/share/perl5:/share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/Perl_Lib/lib64/perl5:/share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/Perl_Lib/CIMS/czplib:$PERL5LIB


cd /share_bio/unisvx1/yangyg_group/chenysh/ZhouQ_TESIS-m6A/NEW_ANALYSIS/P12-WT-miCLIP-Total/bwa
mkdir CTK_Procedure
cd CTK_Procedure


#####1. Parsing SAM file
perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file P12-WT-miCLIP.mutation.txt ../P12-WT-miCLIP.sam P12-WT-miCLIP.tag.bed 


####2. Collapsing PCR duplicates
perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/tag2collapse.pl -v -big --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name P12-WT-miCLIP.tag.bed P12-WT-miCLIP.tag.uniq.bed 
perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/selectRow.pl -q 3 -f 3 P12-WT-miCLIP.mutation.txt P12-WT-miCLIP.tag.uniq.bed > P12-WT-miCLIP.tag.uniq.mutation.txt 


###3. After getting the unique tags of each library, one might concatenate biological replicates, which are distinguished by different colors
###3. 作图用,将一个样本的tag用一种颜色标记
perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/bed2rgb.pl -v -col "128,0,0" P12-WT-miCLIP.tag.uniq.bed P12-WT-miCLIP.tag.uniq.rgb.bed 


###4. As a diagnostic step, get the length distribution of unique tags, which should be a more faithful representation of the library:
awk '{print $3-$2}' P12-WT-miCLIP.tag.uniq.bed | sort -n | uniq -c | awk '{print $2"\t"$1}' > P12-WT-miCLIP.uniq.len.dist.txt


awk '{if($9==">") {print $0}}' P12-WT-miCLIP.tag.uniq.mutation.txt | cut -f 1-6 > P12-WT-miCLIP.tag.uniq.sub.bed
awk '{if($9=="-") {print $0}}' P12-WT-miCLIP.tag.uniq.mutation.txt | cut -f 1-6 > P12-WT-miCLIP.tag.uniq.del.bed
awk '{if($9=="+") {print $0}}' P12-WT-miCLIP.tag.uniq.mutation.txt | cut -f 1-6 > P12-WT-miCLIP.tag.uniq.ins.bed




########################################################
########################################################
cd /share_bio/unisvx1/yangyg_group/chenysh/ZhouQ_TESIS-m6A/NEW_ANALYSIS/P12-WT-miCLIP-Total/bwa/CTK_Procedure
mkdir CIMS_CT
cd CIMS_CT

awk '($8=="C" && $10=="T" && $6=="+") || ($8=="G" && $10=="A" && $6=="-")' ../P12-WT-miCLIP.tag.uniq.mutation.txt |cut -f 1-6 > P12-WT-miCLIP.tag.uniq.sub-CT.bed


cd /share_bio/unisvx1/yangyg_group/chenysh/ZhouQ_TESIS-m6A/NEW_ANALYSIS/P12-WT-miCLIP-Total/bwa/CTK_Procedure/CIMS_CT
###Mutation Mode
perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/CIMS.pl -n 10 -p -v --keep-cache -c cache_mut-CT ../P12-WT-miCLIP.tag.uniq.bed P12-WT-miCLIP.tag.uniq.sub-CT.bed P12-WT-miCLIP.tag.uniq.mut-CT.CIMS.txt
sed '1d' P12-WT-miCLIP.tag.uniq.mut-CT.CIMS.txt|awk '$8>=2 && $8/$7>=0.01 && $8/$7<=0.5'|cut -f1-6 > filter_P12-WT-miCLIP.tag.uniq.mut-CT.CIMS.txt

awk -v OFS="\t" '{if($6=="+"){print $1,$2-3,$3+1,$4,$5,$6}else{print $1,$2-1,$3+3,$4,$5,$6}}' filter_P12-WT-miCLIP.tag.uniq.mut-CT.CIMS.txt|fastaFromBed -fi /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Mouse/ENSEMBLE_68/index/genome.fa  -bed - -s -fo ./filter_P12-WT-miCLIP.tag.uniq.mut-CT.fasta

awk '{if(NR%2==1){tmp=$0}else{split($1,a,"");if(a[3]=="A"){print tmp; print $0}}}' filter_P12-WT-miCLIP.tag.uniq.mut-CT.fasta > m6A_P12-WT-miCLIP.tag.uniq.mut-CT.fasta
sed 's/>//' m6A_P12-WT-miCLIP.tag.uniq.mut-CT.fasta|sed 's/:/\t/'|sed 's/-/\t/'|sed 's/(/\t/' |sed 's/)//'|awk -v OFS="\t" '{if(NR%2==1){chr=$1; chrStart=$2; chrEnd=$3; strand=$4;}else{print chr, chrStart, chrEnd, $0, strand}}' > m6A_P12-WT-miCLIP.tag.uniq.mut-CT.bed





########################################################
##CIMS（mutation model）
########################################################
cd /share_bio/unisvx1/yangyg_group/chenysh/ZhouQ_TESIS-m6A/NEW_ANALYSIS/P12-WT-miCLIP-Total/bwa/CTK_Procedure
mkdir CIMS
cd CIMS

perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/CIMS.pl -n 10 -p -v --keep-cache -c cache_mut ../P12-WT-miCLIP.tag.uniq.bed ../P12-WT-miCLIP.tag.uniq.sub.bed P12-WT-miCLIP.tag.uniq.mut.CIMS.txt
awk '{if($9<=0.05) {print $0}}' P12-WT-miCLIP.tag.uniq.mut.CIMS.txt|sort -k 9,9n -k 8,8nr -k 7,7n > P12-WT-miCLIP.tag.uniq.mut.CIMS.p05.txt
cut -f 1-6 P12-WT-miCLIP.tag.uniq.mut.CIMS.p05.txt > P12-WT-miCLIP.tag.uniq.mut.CIMS.p05.bed

sed '1d' P12-WT-miCLIP.tag.uniq.mut.CIMS.p05.txt|awk -v OFS="\t" '$8>=2 && $8/$7>=0.01 && $8/$7<=0.5{if($6=="+"){print $1,$2-4,$3+4,$4,$5,$6}else{print $1,$2-4,$3+4,$4,$5,$6}}'|cut -f1-6|fastaFromBed -fi /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Mouse/ENSEMBLE_68/index/genome.fa -bed - -s -fo P12-WT-m-miCLIP.tag.uniq.mut.fasta

awk -v OFS="\t" 'NR%2==1{tmp=$1}NR%2==0{if($1~/[G,A][G,A]AC[A,T,C]/){a=$1; gsub("[G,A][G,A]AC[A,T,C]", "\t", a); split(a,b,"\t"); split(b[1],c,""); l=length(c); print tmp, l, $1}}' P12-WT-m-miCLIP.tag.uniq.mut.fasta|sed 's/>//'|sed 's/:/\t/'|sed 's/-/\t/'|sed 's/(/\t/'|sed 's/)//'|awk -v OFS="\t" '$4~/+/{print $1,$2+$5, $2+$5+5, "mut_"$1"_"$2+4"_"$3-4"_"$4, "mut_"$6, $4}$4~/-/{print $1, $3-$5-5, $3-$5, "mut_"$1"_"$2+4"_"$3-4"_"$4, "mut_"$6, $4}' > m6A_P12-WT-m-miCLIP.bed



########################################################
##CITS（truncation model））
########################################################

cd /share_bio/unisvx1/yangyg_group/chenysh/ZhouQ_TESIS-m6A/NEW_ANALYSIS/P12-WT-miCLIP-Total/bwa/CTK_Procedure
mkdir CITS
cd CITS

perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/removeRow.pl -q 3 -f 3 -v ../P12-WT-miCLIP.tag.uniq.bed ../P12-WT-miCLIP.tag.uniq.del.bed > P12-WT-miCLIP.tag.uniq.clean.bed
perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/bedExt.pl -n up -l "-1" -r "-1" -v P12-WT-miCLIP.tag.uniq.clean.bed P12-WT-miCLIP.tag.uniq.clean.trunc.bed  

perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/tag2cluster.pl -big -s -maxgap "-1" -of bed -v ../P12-WT-miCLIP.tag.uniq.bed P12-WT-miCLIP.tag.uniq.cluster.0.bed
awk '{if($5>2) {print $0}}' P12-WT-miCLIP.tag.uniq.cluster.0.bed > P12-WT-miCLIP.tag.uniq.cluster.bed 

perl /share_bio/unisvx1/yangyg_group/chenysh/SOFTWARE/ctk-1.0.3/tag2peak.pl -big -ss -v --prefix "CITS" -gap 25 -p 0.05 -gene P12-WT-miCLIP.tag.uniq.cluster.bed P12-WT-miCLIP.tag.uniq.clean.trunc.bed P12-WT-miCLIP.tag.uniq.clean.CITS.p05.bed 



