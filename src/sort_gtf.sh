#!/bin/bash
if [ $# -ne 2 ]; then
    echo "Usage: $0 in.unsort.gtf out.sort.gtf"
    echo "       sort GTF file based on 'trans' lines"
    exit
fi

# 1. cat head lines to 
#awk '{if($0 ~ /^#/) {print} else {exit}}' $1 > $2
# 2. tag lines based on 'gene' lines
awk '
BEGIN {
chr=0; start=0; chr_m=25; OFS="\t";
chrom["chr1"]=1; chrom["chr2"]=2; chrom["chr3"]=3; chrom["chr4"]=4; chrom["chr5"]=5;
chrom["chr6"]=6; chrom["chr7"]=7; chrom["chr8"]=8; chrom["chr9"]=9; chrom["chr10"]=10;
chrom["chr11"]=11; chrom["chr12"]=12; chrom["chr13"]=13; chrom["chr14"]=14; chrom["chr15"]=15;
chrom["chr16"]=16; chrom["chr17"]=17; chrom["chr18"]=18; chrom["chr19"]=19; chrom["chr20"]=20;
chrom["chr21"]=21; chrom["chr22"]=22; chrom["chrX"]=23; chrom["chrY"]=24; chrom["chrM"]=25;
} 
($0 !~ /^#/ && ($3 ~ "transcript" || $3 == "exon")){ 
    if ($3 == "transcript") {
        if (chrom[$1] == "") {
            chr_m += 1;
            chrom[$1] = chr_m;
        }
        chr = chrom[$1]; start = $4; end = $5;
    }
    print chr, start, end, NR, $0
} ' $1 | sort -n -k1 -n -k2 -n -k3 -n -k4 -n -k5 | awk 'BEGIN{FS="\t"; OFS="\t"} {print $5,$6,$7,$8,$9,$10,$11,$12,$13}' > $2
