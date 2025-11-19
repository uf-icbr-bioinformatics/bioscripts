#!/bin/bash

# Given a BAM file, divide each chromosome in consecutive windows of N basepairs
# (default: 100kb) and create a BAM file containing the coverage of each window.
# Execution is parallelized on the number of chromosomes. 
# Requires samtools in PATH.

set -e

WINSIZE=100000
OUT=""

function usage() {
  cat <<EOF
windowcov.sh - report coverage in consecutive genomic windows

Usage: windowcov.sh [options] file.bam

where options are:

  -w W | Use windows of size W (default: ${WINSIZE})
  -o O | Write output to file O (default: basename of BAM file with .cov.bed extension).

EOF
}

while getopts "w:o:h" opt; do
    case $opt in
	h)
	    usage
	    exit 1
	    ;;
	w)
	    WINSIZE=$OPTARG
	    ;;
	o)
	    OUT=$OPTARG
	    ;;
    esac
done
shift $((OPTIND-1))

BAM=$1
if [[ -z $OUT ]];
then
  OUT=${BAM%.bam}.cov.bed
fi

TMPDIR=$(mktemp -d -p .)

samtools idxstats $BAM | cut -f 1,2 | grep -v \* | grep -v _ | grep -v ERCC > ${TMPDIR}/chrom.sizes

while read chrom maxsize; do
  start=1
  end=${WINSIZE}
  echo $chrom
  while true; do
    echo -e "${chrom}\t${start}\t${end}"
    start=$((start+WINSIZE))
    end=$((end+WINSIZE))
    if [[ $end -gt $maxsize ]]; then
	break
    fi
  done > ${TMPDIR}/${chrom}.bed
  out=${TMPDIR}/${chrom}.cov.bed
  samtools bedcov ${TMPDIR}/${chrom}.bed $BAM > $out &
done < ${TMPDIR}/chrom.sizes
wait

rm -f $OUT
while read chrom maxsize; do
  cat ${TMPDIR}/${chrom}.cov.bed >> $OUT
done < ${TMPDIR}/chrom.sizes

rm -fr ${TMPDIR}
