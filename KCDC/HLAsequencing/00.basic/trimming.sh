#!/bin/bash
#Usage: ./trimming.sh [-i <input_file_list>] [-o <output_prefix>] [-t] [-p]
# Creating fastq file list
#/BDATA/smkim/HLA_seq/downsampling/10X/1003_S71_L002_R1_001_paired.10X.fastq.gz,/BDATA/smkim/HLA_seq/downsampling/10X/1003_S71_L002_R2_001_paired.10X.fastq.gz
#/BDATA/smkim/HLA_seq/downsampling/10X/1006_S72_L002_R1_001_paired.10X.fastq.gz,/BDATA/smkim/HLA_seq/downsampling/10X/1006_S72_L002_R2_001_paired.10X.fastq.gz
#join -1 1 -2 1 <(ls *R1*.fastq.gz | awk '{split($1,arr,"_"); print arr[1],$1}' | sort) <(ls *R2*.fastq.gz | awk '{split($1,arr,"_"); print arr[1],$1}' | sort) | awk -v pwd="$(pwd)" '{print pwd"/"$2","pwd"/"$3}' > fastqs.list
#mkdir trimmed
usage() { echo "Usage: $0 [-i <input_file_list>] [-o <output_prefix>] [-t] [-p]" 1>&2; exit 1; }

while getopts ":i:o:t:p:h" arg; do
	case $arg in
		i)
			i=${OPTARG}
			;;
		o)
			o=${OPTARG}
			;;
		t)
			t=${OPTARG}
			;;
		p)
			p=${OPTARG}
			;;
		h | *)
			usage
			;;
	esac
done
shift $((OPTIND-1))

if [ -z "${i}" ] || [ -z "${o}" ]; then
	usage
fi

echo "i = ${i}"
echo "o = ${o}"
if [ -z "${t}" ]; then
	t=1
fi
echo "t = ${t}"
if [ -z "${p}" ]; then
	p=1
fi
echo "p = ${p}"

TRIMMOMATIC_HOME=/ADATA/dongmun/tools/Trimmomatic-0.39
TRIMMOMATIC=$TRIMMOMATIC_HOME/trimmomatic-0.39.jar
FASTQ_LIST=$i
OUT_DIR=${o%%/}
ADAPTER_DIR=$TRIMMOMATIC_HOME/adapters
ADAPTER=$ADAPTER_DIR/TruSeq3-PE-2.fa:2:30:10

args=

while read line
do
	echo $line

	IFS=',' read -a array <<< $line

	INPUT1=${array[0]}
	INPUT2=${array[1]}

	PAIRED_OUTPUT1=
	UNPAIRED_OUTPUT1=
	PAIRED_OUTPUT2=
	UNPAIRED_OUTPUT2=

	# Output path
	if [[ "$INPUT1" == *.fq.gz ]];then
		PAIRED_OUTPUT1=$OUT_DIR/`basename ${INPUT1%%.*} fq.gz`_paired.fq.gz
		UNPAIRED_OUTPUT1=$OUT_DIR/`basename ${INPUT1%.*} fq.gz`_unpaired.fq.gz
	elif [[ "$INPUT1" == *.fastq.gz ]];then
		PAIRED_OUTPUT1=$OUT_DIR/`basename ${INPUT1%%.*} fastq.gz`_paired.fastq.gz
		UNPAIRED_OUTPUT1=$OUT_DIR/`basename ${INPUT1%%.*} fastq.gz`_unpaired.fastq.gz
	else
		# echo "${INPUT1} file extension is not valid"
		exit
	fi

	if [[ "$INPUT2" == *.fq.gz ]];then
		PAIRED_OUTPUT2=$OUT_DIR/`basename ${INPUT2%%.*} fq.gz`_paired.fq.gz
		UNPAIRED_OUTPUT2=$OUT_DIR/`basename ${INPUT2%%.*} fq.gz`_unpaired.fq.gz
	elif [[ "$INPUT2" == *.fastq.gz ]];then
		PAIRED_OUTPUT2=$OUT_DIR/`basename ${INPUT2%%.*} fastq.gz`_paired.fastq.gz
		UNPAIRED_OUTPUT2=$OUT_DIR/`basename ${INPUT2%%.*} fastq.gz`_unpaired.fastq.gz
	else
		# echo "${INPUT2} file extension is not valid"
		exit
	fi

	# Check if output exists
	if [ -f $PAIRED_OUTPUT1 ]; then
		# echo "'${PAIRED_OUTPUT1}' exists."
		continue
	fi

	if [ -f $PAIRED_OUTPUT2 ]; then
		# echo "'${PAIRED_OUTPUT2}' exists."
		continue
	fi

	arg="${INPUT1} ${INPUT2} ${PAIRED_OUTPUT1} ${UNPAIRED_OUTPUT1} ${PAIRED_OUTPUT2} ${UNPAIRED_OUTPUT2} ILLUMINACLIP:${ADAPTER} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
	args+=$arg'\n'
	#echo $comm
	#args+=($arg)
	#echo "${args[@]}"
done < $FASTQ_LIST

echo -e "$args" | xargs -I{} -n1 -P"$p" sh -c "java -jar ${TRIMMOMATIC} PE -threads ${t} -phred33 {}" >> "trimmomatic-log.txt" 2>&1
