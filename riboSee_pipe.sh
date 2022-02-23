#!/bin/bash
usage () {
  echo "
  This pipeline is soley for bacterial riboSeq (it doesnt account for introns or alt. splicing, otherwise use hisat or star aligners).
  Assumes single end reads (or one paired read).
  Sequences are first aligned to rRNAs and tRNAs annotated in the GTF, and then the unaligned reads are aligned to the reference genome.
  
  Mode 1:
    -dl|--container = download the singularity container to this path and exit
  
  Mode 2:
    -mt|--get_metrics = supply a dir and get metrics for all analyses in that dir, all other options will be ignored if this is non-empyt
  
  Mode 3: 
    -t|--threads
    -g|--genome_reference = path to genome reference
    -gtf|--GTF_reference
    -rd|--read_dir
    -r1|--read2 = the fastq file
    -o|--out_dir
    -n|--nme output name
    -s|--strand = stranded library (yes|no|reverse)
    -sd|--script_directory
    -fq|--fastQC = run fastqc?
    -tm|--trim_fasta = path to multifasta with adapters, linkers etc to trim
    -mm|--max_missmatch (defult = 2)
    -mn|--min_len (defult = 24)
    -mx|--max_len (defult = 36)
  "
}



declare_globals () {
    # if same thing twice will take second one
    while [[ "$#" -gt 0 ]]
    do
        case $1 in
        -t|--threads)
        threads="$2"
        ;;
        -g|--genome_reference)
        ref="$2"
        ;;
        -gtf|--GTF_reference)
        gtf="$2"
        ;;
        -rd|--read_dir)
        read_dir="$2"
        ;;
        -r|--reads) 
        reads="$2"
        ;;
        -o|--out_dir) 
        out_dir="$2"
        ;;
        -n|--name) 
        nme="$2"
        ;;
        -s|--strand) #stranded library (yes|no|reverse)
        strand="$2"
        ;;
        -sd|--script_directory)
        Script_dir="$2"
        ;;
        -dl|--container)
        container="$2"
        ;;
        -mt|--get_metrics)
        get_metrics="$2"
        ;;
        -fq|--fastQC)
        fastQC="Y"
        ;;
        -mm|--max_missmatch)
        max_missmatch="$2"
        ;;
        -mn|--min_len)
        min_len="$2"
        ;;
        -mx|--max_len)
        max_len="$2"
        ;;
        -tm|--trim_fasta)
        trim_fasta="$2"
        ;;
    esac
        shift
    done
}



if [ $# == 0 ]
then
    usage
    exit 1
fi


###############
## Functions ##
qc_trim_SE () {
    java -jar "$TRIM" SE \
      -threads $1 \
      "$2" \
      "${2/.f*/.trimmed.fq.gz}" \
      ILLUMINACLIP:"$3":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:3:10 MINLEN:$4

    #FastQC post
    if  [[ ! -z $fastQC ]]
    then
      fastqc -t $1 "${2/.f*/.trimmed.fq.gz}" -o ./
    fi

    reads="${2/.f*/.trimmed.fq.gz}"
    export reads
}

###############

#####################
## Setup variables ##
declare_globals "$@"

TRIM=/usr/bin/Trimmomatic-0.39/trimmomatic-0.39.jar

if [[ ! -z $container ]]
then
    cd $container
    singularity pull library://semiquant/default/riri
    exit 0
fi

if [[ ! -z $get_metrics ]]
then
    cd $get_metrics
    multiqc "$get_metrics" -n $(basename $get_metrics)
    exit 0
fi

# Set defults
Script_dir_tmp="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
Script_dir="${Script_dir:-$Script_dir_tmp}"
ram=$(expr $threads \* 2)
jav_ram=$(echo "scale=2; $ram*0.8" | bc)
export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
strand="${strand:-no}"

max_missmatch="${max_missmatch:-2}"
min_len="${strand:-24}"
max_len="${strand:-36}"

trim_tmp="${Script_dir}/references/adapts.fasta"
trim_fasta="${trim_fasta:-trim_tmp}"
ref_tmp="${Script_dir}/references/Mycobacterium_tuberculosis_H37Rv_genome_v4.fasta"
ref="${ref:-ref_tmp}"
gtf_tmp="${Script_dir}/references/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff"
gtf="${gtf:-gtf_tmp}"


# Check if programs installed
command -v fastqc >/dev/null 2>&1 || { echo >&2 "I require fastqc but it's not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed. Aborting."; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools but it's not installed. Aborting."; exit 1; }
command -v htseq-count >/dev/null 2>&1 || { echo >&2 "I require htseq but it's not installed. Aborting."; exit 1; }
command -v bowtie >/dev/null 2>&1 || { echo >&2 "I require bowtie but it's not installed. Aborting."; exit 1; }
command -v metagene >/dev/null 2>&1 || { echo >&2 "I require plastid but it's not installed. Aborting."; exit 1; }


# easiest, but not robust
mkdir -p "${read_dir}/${nme}_riri"
cd "${read_dir}/${nme}_riri"
ln -s "$reads" "${read_dir}/${nme}_riri/"
#####################



######################
## setup references ##
if [ ! -e "${1/.f*/}.ebwtl" ] 
then
  bowtie-build --threads $threads "$ref" "${ref/.f*/}"
fi

# generate metagene `roi` file
if [ ! -e ?? ] 
then
  metagene generate "${ref/.f*/}" \
    --landmark cds_start \
    --annotation_files "$gtf"

fi

# split rRNAs an tRNAs
awk '{if ($3 == "tRNA" || $3 == "rRNA") print $0;}' "$gtf" > "${gtf/.g*/_rRNA.gtf}"
ref_rem="${ref/.f*/.rRNAs.fasta}"
bedtools getfasta -fi "$ref" -bed "${gtf/.g*/_rRNA.gtf}" > "$ref_rem"
######################


##############
## Pipeline ##
# Fastqc
if [[ ! -e "${reads/.f*/_fastqc.zip}" ]] && [[ ! -z $fastQC ]]
then
  fastqc -t $threads "$reads" -o "./"
fi

# trim reads
qc_trim_SE $threads "$reads" "$trim_fasta" "$min_len"

# align reads
bowtie --threads $threads --seed 1987 -x "$ref_rem" -q "$reads" -a -v $max_missmatch \
  --un "${reads}_cleaned.fq" > /dev/null

reads="${reads}_cleaned.fq"
bowtie --threads $threads -S --seed 1987 -x "$ref" -q "$reads" -a -v $max_missmatch "${nme}.sam"


# count alignments
htseq-count --nprocesses $threads --type "gene" --idattr "Name" --order "name" --mode "union" --stranded "$strand" \
  --minaqual 10 --nonunique none -f bam "${nme}.bam" "$gtf" --counts_output "${nme/.bam/_HTSeq.counts.tsv}"
# featureCounts -F "GTF" -d 30 -s "$stran_fc" -t "gene" -g "Name" -O -Q 5 --ignoreDup -T $5 -a "$4" "$fCount" -o "${3/.bam/.featCount.counts}" "$3"

# 5′ mapped sites of RPFs
psite "${ref/.f*/.txt}" "${nme}_riboprofile" \
  --min_length $min_len \
  --max_length $max_len \
  --require_upstream \
  --count_files "${nme}.bam"

##############        









