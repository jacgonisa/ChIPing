#!/bin/bash

PARAMS=$1
if [ $# -ne 1 ]
then
	echo "Number of arguments is: $#"
 	echo "Please, only 1 argument for this function!"
	echo "Usage chipipe.sh <params.file>"
	exit
fi

PARAMS=$1

INSDIR=$(grep installation $PARAMS | awk '{print($2)}')
echo "Installation directory is: $INSDIR"

WD=$(grep working $PARAMS | awk '{print($2)}')
echo "Working directory is: $WD"

EXP=$(grep experiment $PARAMS | awk '{print($2)}')
echo "Experiment directory is: $EXP"

NUMREPLICAS=$(grep number_replicas $PARAMS | awk '{print($2)}')
echo " Number of replica: $NUMREPLICAS"

GENOME=$(grep genome $PARAMS | awk '{print($2)}')
echo "Genome directory is: $GENOME"

ANNOTATION=$(grep annotation $PARAMS | awk '{print($2)}')
echo "Annotation directory is: $ANNOTATION"

CHR=$(grep chromosomes $PARAMS | awk '{print($2)}')
echo "Universe of chromosomes directory is: $CHR"

PEAK=$(grep peak $PARAMS | awk '{print($2)}')
echo "Peak type: $PEAK"

SINGLE=$(grep single $PARAMS | awk '{print($2)}')
echo "SINGLE or PAIRED: $SINGLE"

TSSUP=$(grep upstream $PARAMS | awk '{print($2)}')
echo "TSS region upstream is: $TSSUP"

TSSDOWN=$(grep downstream $PARAMS | awk '{print($2)}')
echo "TSS region downstream is: $TSSDOWN"

CHIPS=()
INPUTS=()
i=0
if [ $SINGLE -eq 1 ]
then
	while [ $i -lt $NUMREPLICAS ]
	do
		j=$(($i + 1))
		CHIPS[$i]=$(grep path_sample_chip_$j $PARAMS | awk '{print($2)}')
		INPUTS[$i]=$(grep path_sample_input_$j $PARAMS | awk '{print($2)}')
		((i++))
	done
elif  [ $SINGLE -eq 2] 
then
	echo "no hay paired todavia"
else
	echo "No allowed input for single/paired end reads determination"
fi

echo "Samples = "
echo "${CHIPS[@]}"
echo "${INPUTS[@]}"

#Generationg work space

echo "====================="
echo "GENERATING WORK SPACE"
echo "====================="

cd $WD
mkdir $EXP
cd $EXP
mkdir genome annotation results samples scripts
cp $GENOME genome/genome.fa
cp $ANNOTATION annotation/annotation.gtf
cd samples

if [ $SINGLE -eq 1 ]
then
	i=1
	while [ $i -le $NUMREPLICAS ]
	do
		mkdir replica_$i
		cd replica_$i
		mkdir chip input replica_results
		j=$(($i-1))
		cp ${CHIPS[$j]} chip/sample_chip_$i.fq.gz
		cp ${INPUTS[$j]} input/sample_input_$i.fq.gz
		cd ..
		((i++))
	done
fi

echo "===================="
echo "WORKSPACE CREATED"
echo "===================="

echo "====================="
echo "Creating index"
echo "====================="

cd ../genome
bowtie2-build genome.fa index
echo "Files size:" du -h *

echo "====================="
echo "Processing samples"
echo "====================="

cd ../results


i=1
while [ $i -le $NUMREPLICAS ]
do
	sbatch --job-name=proc_$i --output=sample_proc_$i --error=err_$i $INSDIR/sample_proc.sh $WD $i $PEAK $NUMREPLICAS $INSDIR $EXP $CHR $TSSUP $TSSDOWN $GENOME
	((i++))
done
