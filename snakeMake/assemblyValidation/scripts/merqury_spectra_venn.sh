#!/bin/bash

function usage {
  echo "Usage: spectra-cn.sh <read.meryl> <asm1.fasta> [asm2.fasta] out-prefix"
  echo -e "\t<read.meryl>\t: Generated with meryl count from i.e. illumina wgs reads"
  echo -e "\t<asm1.fasta>\t: haplotype 1 assembly. gzipped or not"
  echo -e "\t<ASM prefix>\t: A prefix to ensure that the temp files are unique"
  echo -e "\t<out-prefix>: output prefix. Required."
}

source $MERQURY/util/util.sh

mkdir -p merqury/$asm

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
read_hist=merqury/$asm/$asm.$db.hist
hist_ploidy=merqury/$asm/$asm.$db.hist.ploidy
read_filt=merqury/$asm/$asm.$db.filt
hist=merqury/$asm/$name.$asm.spectra-cn.hist
hist_asm_only=merqury/$asm/$name.$asm.only.hist

asm1=$asm
asm_db=${asm}.meryl
read_k_copy0=merqury/$asm/$name.read.k$k.$asm.0.meryl
read_k_copy1=merqury/$asm/$name.read.k$k.$asm.1.meryl
read_k_copy2=merqury/$asm/$name.read.k$k.$asm.2.meryl
read_k_copy3=merqury/$asm/$name.read.k$k.$asm.3.meryl
read_k_copy4=merqury/$asm/$name.read.k$k.$asm.4.meryl
read_k_copygt4=merqury/$asm/$name.read.k$k.$asm.gt4.meryl
read_asm_int=merqury/$asm/$name.read.k$k.$asm1.meryl

asm_only=merqury/$asm/$name.$asm.0.meryl
asm_solid=merqury/$asm/$asm.solid.meryl

# outputs
asm_qv=merqury/$asm/$name.qv
perseq_qv=merqury/$asm/$name.$asm.qv
completeness=merqury/$asm/$name.completeness.stats
hist=merqury/$asm/$name.spectra-asm.hist
hist_asm_dist_only=merqury/$asm/$name.dist_only.hist

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
      meryl intersect output merqury/$asm/$name.read.k$k.$asm.$i.meryl $read [ equal-to $i $asm_db ]
      meryl histogram merqury/$asm/$name.read.k$k.$asm.$i.meryl | awk -v cn=$i '{print cn"\t"$0}' >> $hist
      #rm -r $name.read.k$k.$asm.$i.meryl
      echo
  done

  echo "Copy >4 .."
  meryl intersect output $read_k_copygt4 $read [ greater-than $i $asm_db ]
  meryl histogram $read_k_copygt4 | awk -v cn=">$i" '{print cn"\t"$0}' >> $hist
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
Rscript $MERQURY/plot/plot_spectra_cn.R -f $hist -o $name.$asm.spectra-cn -z $hist_asm_only --pdf"
Rscript $MERQURY/plot/plot_spectra_cn.R -f $hist -o $name.$asm.spectra-cn -z $hist_asm_only --pdf
echo

echo "# QV statistics"
ASM_ONLY=`meryl statistics $asm_only  | head -n4 | tail -n1 | awk '{print $2}'`
      TOTAL=`meryl statistics $asm_db  | head -n4 | tail -n1 | awk '{print $2}'`
      ERROR=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (1-(1-$1/$2)^(1/k))}'`
      QV=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'`
      echo -e "$asm\t$ASM_ONLY\t$TOTAL\t$QV\t$ERROR" >> $asm_qv
echo

echo "# Per seq QV statistics"
meryl-lookup -existence -sequence $asm_fa -mers $asm_only | \
awk -v k=$k '{print $1"\t"$NF"\t"$(NF-2)"\t"(-10*log(1-(1-$NF/$(NF-2))^(1/k))/log(10))"\t"(1-(1-$NF/$(NF-2))^(1/k))}' > $perseq_qv
echo

      echo "# k-mer completeness (recovery rate) with solid k-mers for $asm with > $filt counts"
meryl intersect output $asm_solid $asm_db $read_solid
      TOTAL=`meryl statistics $read_solid | head -n3 | tail -n1 | awk '{print $2}'`
      ASM=`meryl statistics $asm_solid | head -n3 | tail -n1 | awk '{print $2}'`
      echo -e "${asm}\tall\t${ASM}\t${TOTAL}" | awk '{print $0"\t"((100*$3)/$4)}' >> $completeness
rm -r $asm_solid
echo

echo "# $asm1 only"
meryl intersect output $read_asm_int $read $asm_db

echo "# Write output"
echo -e "Assembly\tkmer_multiplicity\tCount" > $hist
meryl histogram $read_k_copy0 | awk '{print "read-only\t"$0}' >> $hist
meryl histogram $read_asm_int | awk -v hap="${asm1}" '{print hap"\t"$0}' >> $hist

echo "# Get asm only for spectra-asm"
ASM1_ONLY=`meryl statistics $asm_only | head -n3 | tail -n1 | awk '{print $2}'`
echo -e "${asm1}\t0\t$ASM1_ONLY" > $hist_asm_dist_only

#echo "#	Plot $hist"
#echo "\
#Rscript $MERQURY/plot/plot_spectra_cn.R -f $hist -o $name.$asm.spectra-asm -z $hist_asm_dist_only --pdf"
#echo

echo "# Cleaning up"
rm -r $read_filt $hist_asm_only $read_k_copy0 $read_k_copy1 $read_k_copy2 $read_k_copy3
rm -r $read_k_copy4 $read_k_copygt4 $asm_only
