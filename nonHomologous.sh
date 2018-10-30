#! /bin/bash

set -e
set -u
set -o pipefail

usage() {
    echo
    echo "Usage: $(basename "$0") [-h] -b BAMDIR -o OUTDIR -- run pipeline to find non-homologous recombination in a set of bam files from a SNP pipeline output." >&2
    echo
    echo "    -h                        Print this help text and exit"
    echo "    -b BAMDIR                 Directory containing the bam files to analyze"
    echo "    -o OUTDIR                 Directory to put the output files in"
    echo
}

while getopts ":b:o:h" opt; do
    case $opt in
        h)
          usage
          exit 0
          ;;
        b)
          bamdir=$OPTARG
          ;;
        o)
          outdir=$OPTARG
          ;;
        \?)
          echo "Invalid option: -$OPTARG" >&2
          usage
          exit 1
          ;;
        :)
          echo "Option -$OPTARG requires an argument." >&2
          usage
          exit 1
          ;;
    esac
done

mkdir -p $outdir
cd $bamdir
assembleUnmappedReads

for fasta in *.fasta; do
    awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' $fasta >> ${outdir}/accessory_genomes_combined.fasta
done

alignAccessoryGenomes -f ${outdir}/accessory_genomes_combined.fasta -l 5000 -d 2 -o $outdir

cd $outdir
for ref in `grep -o -E "^>\S+" accessory_genomes_trimmed.fasta | tr -d ">"`; do 
    for query in `grep -o -E "^>\S+" accessory_genomes_trimmed.fasta | tr -d ">"`; do 
        show-aligns -r -w 80 pairwise_filtered.delta $ref $query >> pairwise_aligned.out 2> /dev/null
    done
done 


