#!/bin/bash

echo "Usage: spectra-cn.sh <read.meryl> <asm1.fasta> [asm2.fasta] out-prefix"
echo -e "\t<read.meryl>\t: Generated with meryl count from i.e. illumina wgs reads"
echo -e "\t<asm1.fasta>\t: haplotype 1 assembly. gzipped or not"
echo -e "\t<ASM prefix>\t: A prefix to ensure that the temp files are unique"
echo -e "\t<out-prefix>: output prefix. Required."
echo -e "\t\tWhen only <asm1.fasta> is given, results will be generated in haploid mode."
echo -e "\t\tWhen <asm2.fasta> is given, results will be generated for each asm1 asm2 haploid assembly and asm1+asm2 diploid assembly."
echo


if [[ $# -lt 3 ]]; then
    echo "No args provided. Exit."
    exit -1
fi

source $MERQURY/util/util.sh

read=`link $1`
asm1_fa=`link $2`
asm=$3
name=$4

k=`meryl print $read | head -n 2 | tail -n 1 | awk '{print length($1)}'`
echo "Detected k-mer size $k"
echo

asm2_fa=""


if [ -s $name ]; then
        echo "$name already exists. Provide a different name."
        exit -1
fi


asm1=$asm

echo "# Get solid k-mers"
### Taken from filt.sh to avoid file collisions
db=$read
db=${db/.meryl}

echo "Generate $db.hist"
meryl histogram $db.meryl > $asm.$db.hist

echo "
java -jar -Xmx1g $MERQURY/eval/kmerHistToPloidyDepth.jar $asm.$db.hist
"
java -jar -Xmx1g $MERQURY/eval/kmerHistToPloidyDepth.jar $asm.$db.hist > $asm.$db.hist.ploidy

cat $asm.$db.hist.ploidy

filt=`sed -n 2p $asm.$db.hist.ploidy | awk '{print $NF}'`

echo "
Filter out kmers <= $filt"

echo "
meryl greater-than $filt output $asm.$db.gt$filt.meryl $db.meryl
"
meryl greater-than $filt output $asm.$db.gt$filt.meryl $db.meryl
echo $filt > $asm.$db.filt

###

#filt=`cat ${read/.meryl/.filt}`
read_solid=$asm.$db.gt$filt.meryl

echo "=== Generate spectra-cn plots per assemblies and get QV, k-mer completeness ==="
echo
asm_fa=$asm1_fa

	#asm=`echo $asm_fa | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

	if [ ! -e $asm.meryl ]; then
	    echo "# Generate meryl db for $asm"
	    meryl count k=$k output ${asm}.meryl $asm_fa
	    echo
	fi

	echo "# Collect read counts per asm copies"
	hist=$name.$asm.spectra-cn.hist
	hist_asm_only=$name.$asm.only.hist

	if [[ -s $hist ]]; then
		echo
		echo "*** $hist found. ***"
		echo
	else

		echo -e "Copies\tkmer_multiplicity\tCount" > $hist

		echo "# Read only"
		meryl difference output $name.read.k$k.$asm.0.meryl $read $asm.meryl
		meryl histogram $name.read.k$k.$asm.0.meryl | awk '{print "read-only\t"$0}' >> $hist

		echo "# Copy 1 ~ 4"
		for i in $(seq 1 4)
		do
		    echo "Copy = $i .."
		    meryl intersect output $name.read.k$k.$asm.$i.meryl $read [ equal-to $i ${asm}.meryl ]
		    meryl histogram $name.read.k$k.$asm.$i.meryl | awk -v cn=$i '{print cn"\t"$0}' >> $hist
		    #rm -r $name.read.k$k.$asm.$i.meryl
		    echo
		done

		echo "Copy >4 .."
		meryl intersect output $name.read.k$k.$asm.gt$i.meryl $read [ greater-than $i ${asm}.meryl ]
		meryl histogram $name.read.k$k.$asm.gt$i.meryl | awk -v cn=">$i" '{print cn"\t"$0}' >> $hist
		#rm -r $name.read.k$k.$asm.gt$i.meryl
		echo
	fi

	echo "# Copy numbers in k-mers found only in asm"
	meryl difference output $name.$asm.0.meryl ${asm}.meryl $read
	PRESENT=`meryl statistics $name.${asm}.0.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
	DISTINCT=`meryl statistics $name.${asm}.0.meryl  | head -n3 | tail -n1 | awk '{print $2}'`
	MULTI=$(($PRESENT-$DISTINCT))
	echo -e "1\t0\t$DISTINCT" > $hist_asm_only
	echo -e "2\t0\t$MULTI" >> $hist_asm_only
	echo

	echo "# Plot $hist"
	echo "\
	Rscript $MERQURY/plot/plot_spectra_cn.R -f $hist -o $name.$asm.spectra-cn -z $hist_asm_only"
	Rscript $MERQURY/plot/plot_spectra_cn.R -f $hist -o $name.$asm.spectra-cn -z $hist_asm_only
	echo

	echo "# QV statistics"
	ASM_ONLY=`meryl statistics $name.${asm}.0.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
        TOTAL=`meryl statistics ${asm}.meryl  | head -n4 | tail -n1 | awk '{print $2}'`
        ERROR=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (1-(1-$1/$2)^(1/k))}'`
        QV=`echo "$ASM_ONLY $TOTAL" | awk -v k=$k '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'`
        echo -e "$asm\t$ASM_ONLY\t$TOTAL\t$QV\t$ERROR" >> merqury/$asm/$name.qv
	echo

	echo "# Per seq QV statistics"
	meryl-lookup -existence -sequence $asm_fa -mers $name.$asm.0.meryl/ | \
	awk -v k=$k '{print $1"\t"$NF"\t"$(NF-2)"\t"(-10*log(1-(1-$NF/$(NF-2))^(1/k))/log(10))"\t"(1-(1-$NF/$(NF-2))^(1/k))}' > $name.$asm.qv
	echo

        echo "# k-mer completeness (recoveray rate) with solid k-mers for $asm with > $filt counts"
	meryl intersect output $asm.solid.meryl $asm.meryl $read_solid
        TOTAL=`meryl statistics $read_solid | head -n3 | tail -n1 | awk '{print $2}'`
        ASM=`meryl statistics $asm.solid.meryl | head -n3 | tail -n1 | awk '{print $2}'`
        echo -e "${asm}\tall\t${ASM}\t${TOTAL}" | awk '{print $0"\t"((100*$3)/$4)}' >> merqury/$asm/$name.completeness.stats
	rm -r $asm.solid.meryl
	echo

	echo "# Generate ${asm}_only.tdf"
	if [[ ! -e "$asm_fa.fai" ]]; then
		echo "# Index $asm_fa"
		samtools faidx $asm_fa
		echo
	fi

	if [[ ! -e "${asm}_only.tdf" ]]; then
		meryl-lookup -dump -memory 6 -sequence $asm_fa -mers $name.${asm}.0.meryl | awk '$(NF-4)=="T" {print $1"\t"$(NF-5)"\t"($(NF-5)+21)}' > ${asm}_only.bed
		igvtools count ${asm}_only.bed ${asm}_only.tdf $asm_fa.fai
		echo "${asm}_only.tdf generated."
	else
		echo
		echo "*** ${asm}_only.tdf found. ***"
	fi
	echo


hist_asm_dist_only=merqury/$asm/$name.dist_only.hist
if [[ "$asm2_fa" = "" ]]; then
	echo "No asm2_fa given. Done."

	hist=merqury/$asm/$name.spectra-asm.hist

	if [[ -s $hist ]]; then
		echo "*** Found $hist ***"
	else
		echo "# $asm1 only"
		meryl intersect output $name.read.k$k.$asm1.meryl $read ${asm1}.meryl

		echo "# Write output"
		echo -e "Assembly\tkmer_multiplicity\tCount" > $hist
		meryl histogram $name.read.k$k.$asm1.0.meryl | awk '{print "read-only\t"$0}' >> $hist
		meryl histogram $name.read.k$k.$asm1.meryl | awk -v hap="${asm1}" '{print hap"\t"$0}' >> $hist

		echo "# Get asm only for spectra-asm"
		ASM1_ONLY=`meryl statistics $name.${asm1}.0.meryl | head -n3 | tail -n1 | awk '{print $2}'`
		echo -e "${asm1}\t0\t$ASM1_ONLY" > $hist_asm_dist_only
	fi

	echo "#	Plot $hist"
	echo "\
	Rscript $MERQURY/plot/plot_spectra_cn.R -f $hist -o $name.spectra-asm -z $hist_asm_dist_only --pdf"
	Rscript $MERQURY/plot/plot_spectra_cn.R -f $hist -o $name.spectra-asm -z $hist_asm_dist_only --pdf
	echo

	echo "# Clean up"
	rm -r ${asm1}.0.meryl read.k$k.$asm1.0.meryl read.k$k.$asm1.meryl $read_solid
	echo "Done!"

	exit 0
fi
