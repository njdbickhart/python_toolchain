#!/usr/bin/env bash
# $1 = depth file
# $2 = vcf file
# $3 = output qv estimate

NUM_BP=`perl -lane 'print $F[0];' < $1`
NUM_SNP=`cat $2 |grep -v "#" | awk -F "\t" '{if (!match($NF, "0/1")) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8}' | tr ';' ' ' | sed s/AB=//g | awk -v WEIGHT=0 '{if ($6 >= WEIGHT) print $0}' | awk -v SUM=0 '{if (length($4) == length($5)) { SUM+=length($4); } else if (length($4) < length($5)) { SUM+=length($5)-length($4); } else { SUM+=length($4)-length($5)}} END { print SUM}'`
echo "num snp: "$NUM_SNP
echo "num bp: "$NUM_BP

perl -e 'chomp(@ARGV); $ns = $ARGV[0]; $nb = $ARGV[1]; print (-10 * log($ns/$nb)/log(10)); print "\n";' $NUM_SNP $NUM_BP > $3
