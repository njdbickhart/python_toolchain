#!/usr/bin/env bash
# $1 = directory of fasta files
# $2 = extension of fasta files
# $3 = output file

for i in `ls $1/*.$2`
do
  echo $i
done > /dev/stderr

for i in `ls $1/*.$2`
do
  cat $i
done > $3
