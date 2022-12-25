#!/bin/bash

#SBATCH --export=ALL

WD=$1
EXP=$2
PEAK=$3
NUMREPLICAS=$4
INSDIR=$5
CHR=$6
TSSUP=$7
TSSDOWN=$8
GENOME=$9

echo "$NUMREPLICAS"
echo "$PEAK"
echo "$EXP"
echo "$GENOME"

cd $WD/$EXP

if [ $PEAK -eq 1 ]
then
	EXT=$(echo "narrowPeak")
elif [ $PEAK -eq 2 ]
then
	EXT=$(echo "broadPeak")
fi

i=3
if [ $NUMREPLICAS -eq 1 ]
then
	mv samples/replica_1/replica_results/1_peaks.${EXT} results/merged_2.${EXT}
else
	bedtools intersect -a samples/replica_1/replica_results/1_peaks.${EXT} -b samples/replica_2/replica_results/2_peaks.${EXT} > results/merged_2.${EXT}
	
	if [ $NUMREPLICAS -ge 3 ]
	then

		while [ $i -le $NUMREPLICAS ] 
		do
			j=$(($i-1))
			bedtools intersect -a results/merged_$((j)).${EXT} -b samples/replica_$i/replica_results/$((i))_peaks.${EXT} > results/merged_$((i)).${EXT}
			((i++))
		done
	fi
fi

cd results
i=$((i-1))

mkdir kegg_images

##The end

Rscript $INSDIR/chip_bash.R merged_$((i)).${EXT} $CHR $TSSUP $TSSDOWN $PEAK

#Homer

findMotifsGenome.pl merged_$((i)).${EXT} $GENOME . -len 9  -size 100
