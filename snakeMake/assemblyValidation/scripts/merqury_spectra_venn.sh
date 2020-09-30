#!/bin/bash

function usage {
  echo "Usage: spectra-cn.sh <read.meryl> <asm1.fasta> [asm2.fasta] out-prefix"
  echo -e "\t<read.meryl>\t: Generated with meryl count from i.e. illumina wgs reads"
  echo -e "\t<asm1.fasta>\t: haplotype 1 assembly. gzipped or not"
  echo -e "\t<ASM prefix>\t: A prefix to ensure that the temp files are unique"
  echo -e "\t<out-prefix>: output prefix. Required."
}

source $MERQURY/util/util.sh

read=`link $1`
asm1_fa=`link $2`
asm=$3
name=$4

if [[ $# -lt 3 ]]; then
    usage
    exit -1
fi

k=`meryl print $read | head -n 2 | tail -n 1 | awk '{print length($1)}'`
echo "Detected k-mer size $k"
echo

### Variables
### Taken from filt.sh to avoid file collisions
db=$read
db=${db/.meryl}
read_hist=$asm.$db.hist
hist_ploidy=$asm.$db.hist.ploidy
read_filt=$asm.$db.filt
hist=$name.$asm.spectra-cn.hist
hist_asm_only=$name.$asm.only.hist

asm_db=${asm}.meryl


echo "# Get solid k-mers"
echo "Generate $db.hist"
meryl histogram $db.meryl > $read_hist

echo "
java -jar -Xmx1g $MERQURY/eval/kmerHistToPloidyDepth.jar $read_hist
"
java -jar -Xmx1g $MERQURY/eval/kmerHistToPloidyDepth.jar $read_hist > $hist_ploidy

cat $hist_ploidy

filt=`sed -n 2p $hist_ploidy | awk '{print $NF}'`

echo "
Filter out kmers <= $filt"

read_gtfilt=$asm.$db.gt$filt.meryl
echo "
meryl greater-than $filt output $read_gtfilt $db.meryl
"
meryl greater-than $filt output $read_gtfilt $db.meryl
echo $filt > $read_filt


read_solid=$asm.$db.gt$filt.meryl

echo "=== Generate spectra-cn plots per assemblies and get QV, k-mer completeness ==="
echo
asm_fa=$asm1_fa


if [ ! -e $asm_db ]; then
    echo "# Generate meryl db for $asm"
    meryl count k=$k output $asm_db $asm_fa
    echo
fi

echo "# Collect read counts per asm copies"

if [[ -s $hist ]]; then
  echo
  echo "*** $hist found. ***"
  echo
else

echo -e "Copies\tkmer_multiplicity\tCount" > $hist
