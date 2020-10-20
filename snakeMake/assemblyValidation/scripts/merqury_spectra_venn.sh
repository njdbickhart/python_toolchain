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
read_k_copy0=$name.read.k$k.$asm.0.meryl
read_k_copy1=$name.read.k$k.$asm.1.meryl
read_k_copy2=$name.read.k$k.$asm.2.meryl
read_k_copy3=$name.read.k$k.$asm.3.meryl
read_k_copy4=$name.read.k$k.$asm.4.meryl
read_k_copygt4=$name.read.k$k.$asm.gt4.meryl

asm_only=$name.$asm.0.meryl

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
  echo "# Read only"
  meryl difference output $read_k_copy0 $read $asm_db
  meryl histogram $read_k_copy0 | awk '{print "read-only\t"$0}' >> $hist

  echo "# Copy 1 ~ 4"
  for i in $(seq 1 4)
  do
      echo "Copy = $i .."
      meryl intersect output $name.read.k$k.$asm.$i.meryl $read [ equal-to $i $asm_db ]
      meryl histogram $name.read.k$k.$asm.$i.meryl | awk -v cn=$i '{print cn"\t"$0}' >> $hist
      #rm -r $name.read.k$k.$asm.$i.meryl
      echo
  done

  echo "Copy >4 .."
  meryl intersect output read_k_copygt4 $read [ greater-than $i $asm_db ]
  meryl histogram read_k_copygt4 | awk -v cn=">$i" '{print cn"\t"$0}' >> $hist
  #rm -r $name.read.k$k.$asm.gt$i.meryl
  echo
fi

echo "# Copy numbers in k-mers found only in asm"
meryl difference output $asm_only $asm_db $read
PRESENT=`meryl statistics $asm_only  | head -n4 | tail -n1 | awk '{print $2}'`
DISTINCT=`meryl statistics $asm_only  | head -n3 | tail -n1 | awk '{print $2}'`
MULTI=$(($PRESENT-$DISTINCT))
echo -e "1\t0\t$DISTINCT" > $hist_asm_only
echo -e "2\t0\t$MULTI" >> $hist_asm_only
echo

echo "# Plot $hist"
echo "\
Rscript $MERQURY/plot/plot_spectra_cn.R -f $hist -o $name.$asm.spectra-cn -z $hist_asm_only"
Rscript $MERQURY/plot/plot_spectra_cn.R -f $hist -o $name.$asm.spectra-cn -z $hist_asm_only
echo
