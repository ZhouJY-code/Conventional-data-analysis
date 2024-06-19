
export PATH=/software/biosoft/software/python/python2021/bin:$PATH
export PATH=/software/biosoft/software/samtools-1.9:$PATH

raw="all.fastq.gz"

path="nanopore"

cd $path

mkdir nanoQC

nanoQC  $raw -o nanoQC
NanoStat --fastq  $raw --outdir statreports

NanoPlot --fastq $raw -t 16  --plots hex dot kde -o nanoplot
NanoPlot --summary sequencing_summary.txt --loglength -o summary

/software/biosoft/software/Filtlong/bin/filtlong  --min_mean_q 7 $raw |  gzip > clean.filtlong.fq.gz


/software/biosoft/software/Filtlong/bin/minimap2-2.24_x64-linux/minimap2 --MD -a -k 13 -x splice -N 32 -un T2T-CoV2.mmi clean.filtlong.fq.gz > alignment.sam
samtools sort -@ 8 -o bam -o sorted.bam alignment.sam
samtools index sorted.bam
samtools faidx ref.fq

