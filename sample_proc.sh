#!/bin/bash

#SBATCH --export=ALL

WD=$1
i=$2
PEAK=$3
NUMREPLICAS=$4
INSDIR=$5
EXP=$6
CHR=$7
TSSUP=$8
TSSDOWN=$9
GENOME=${10}

echo "$GENOME"

cd $WD/$EXP/samples/replica_$i

cd chip

if [ -f sample_chip_${i}_2.fq.gz ]
then
	echo "Espero que no pase"
else
	fastqc sample_chip_$i.fq.gz
	bowtie2 -x ../../../genome/index -U sample_chip_$i.fq.gz -S chip_$i.sam 
fi


cd ../input

if [ -f sample_input_${i}_2.fq.gz ]
then
        echo "Espero que no pase"
else
        fastqc sample_input_$i.fq.gz
        bowtie2 -x ../../../genome/index -U sample_input_$i.fq.gz -S input_$i.sam
fi

samtools sort -o input_$i.bam input_$i.sam
rm input_$i.sam
samtools index input_$i.bam

cd ../chip
samtools sort -o chip_$i.bam chip_$i.sam
rm chip_$i.sam
samtools index chip_$i.bam

cd ../replica_results

# Peak calling

if [ $PEAK -eq 1 ]
then
	echo "PEAK_1"
	macs2 callpeak -t ../chip/chip_$i.bam -c ../input/input_$i.bam -f BAM -n $i
elif [ $PEAK -eq 2 ]
then
	macs2 callpeak --broad -t ../chip/chip_$i.bam -c ../input/input_$i.bam -f BAM -n $i
fi

echo "Peak calling $i done!" >> ../../../results/blackboard

cd ../../../

echo "==============="
echo "REPLICA $i DONE"
echo "==============="


NUMPROC=$(wc -l results/blackboard | awk '{print($1)}')

if [ $NUMPROC -eq $NUMREPLICAS ]
then
	sbatch --job-name=sample_proc_$i --output=out_$i --error=err_$i $INSDIR/peak_call.sh $WD $EXP $PEAK $NUMREPLICAS $INSDIR $CHR $TSSUP $TSSDOWN $WD/$EXP/genome/genome.fa
fi


