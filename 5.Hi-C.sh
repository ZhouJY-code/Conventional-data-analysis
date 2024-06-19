
export PATH=/software/biosoft/software/python/python2020/bin:$PATH
source activate gridtools
export PATH=/software/biosoft/software/samtools-1.9:$PATH

path="HiC"

cd $path

i="LQL220417-10"
juicer_path="./HiC/juicer"

$juicer_path/scripts/juicer.sh -g hg38-cov-rDNA -d ${i} -s DpnII -p $juicer_path/references/hg38-cov-rDNA.size.txt -y $juicer_path/restriction_sites/hg38-cov-rDNA_DpnII.txt -z $juicer_path/references/hg38-cov-rDNA.fasta -D $juicer_path -t 30
