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
    -r|--reads = the fastq file
    -o|--out_dir
    -n|--name = output name
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
        -ca|--cut_adapt)
        cut_adapt="$2"
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
    echo "Started trimming"
    out_fq="${2/.f*/.trimmed.fq.gz}"
    out_fq="$(basename $out_fq)"
    out_fq="${5}/${out_fq}"

    if [[ ! -z $cut_adapt ]]
    then
      cutadapt --cores=$1 \
        --quality-cutoff 10,10 \
        --minimum-length $4 \
        -g $cut_adapt \
        -o "${5}/${out_fq}" "$2"
      # Regular 3’ adapter  -a ADAPTER
      # Regular 5’ adapter  -g ADAPTER
    else
      java -jar "$TRIM" SE \
        -threads $1 \
        "$2" \
        "$out_fq" \
        ILLUMINACLIP:"$3":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:3:10 MINLEN:$4
    fi

    #FastQC post
    if  [[ ! -z $fastQC ]]
    then
      echo "FastQC post trimming"
      fastqc -t $1 "${2/.f*/.trimmed.fq.gz}" -o "${5}/"
    fi

    reads="$out_fq"
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
    singularity pull library://semiquant/default/riri:v0.1
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
threads="${threads:-4}"
ram=$(expr $threads \* 2)
jav_ram=$(echo "scale=2; $ram*0.8" | bc)
export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
strand="${strand:-no}"

max_missmatch="${max_missmatch:-2}"
min_len="${min_len:-24}"
max_len="${max_len:-36}"

trim_tmp="${Script_dir}/references/adapts.fasta"
trim_fasta="${trim_fasta:-trim_tmp}"
ref_tmp="${Script_dir}/references/NC_000962.fasta"
ref="${ref:-ref_tmp}"
gtf_tmp="${Script_dir}/references/NC_000962.gff"
gtf="${gtf:-gtf_tmp}"


# Check if programs installed
command -v fastqc >/dev/null 2>&1 || { echo >&2 "I require fastqc but it's not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed. Aborting."; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools but it's not installed. Aborting."; exit 1; }
command -v htseq-count >/dev/null 2>&1 || { echo >&2 "I require htseq but it's not installed. Aborting."; exit 1; }
command -v bowtie >/dev/null 2>&1 || { echo >&2 "I require bowtie but it's not installed. Aborting."; exit 1; }
command -v metagene >/dev/null 2>&1 || { echo >&2 "I require plastid but it's not installed. Aborting."; exit 1; }


# easiest, but not robust
mkdir -p "${out_dir}/${nme}_riri"
cd "${out_dir}/${nme}_riri"
ln -s "$reads" "${out_dir}/${nme}_riri/"
out_dir="${out_dir}/${nme}_riri/"
#####################



######################
## setup references ##
if [ ! -f "${ref}.1.ebwt" ] 
then
  bowtie-build --threads $threads "$ref" "$ref"
fi

ref_rem="${ref/.f*/_rRNAs.fasta}"
if [ ! -f "${ref_rem}.1.ebwt" ] 
then
  # split rRNAs an tRNAs
  awk '{if ($3 == "tRNA" || $3 == "rRNA") print $0;}' "$gtf" > "${gtf/.g*/_rRNA.gtf}"
  bedtools getfasta -fi "$ref" -bed "${gtf/.g*/_rRNA.gtf}" > "$ref_rem"
  bowtie-build --threads $threads "$ref_rem" "$ref_rem"
fi


# generate metagene `roi` file
if [ ! -f "${ref/.f*/_rois.txt}" ]
then
    if [ ! -f "${gtf_ps}_5p.gtf" ]
    then
        if [[ "${gtf##*.}" == "gff" ]]
        then
          if [ -f "${gtf/.gff/.gtf}" ]
          then
            gtf_ps="${gtf/.gff/.gtf}"
          else
            gffread "$gtf" -T -o "${gtf/.gff/_tmp.gtf}"
            gtf_ps="${gtf/.gff/_tmp.gtf}"
            fi
          fi
        # run python script to make understandable gtf_tmp
        # this is only bcause bacterial annotations in file
        python3 "${Script_dir}/gtf_primer.py" --gtf_in "$gtf_ps" > /dev/null 2>&1
    fi
    metagene generate "${ref/.f*/}" \
      --landmark cds_start \
      --annotation_files "${gtf_ps}_5p.gtf"
fi

######################


##############
## Pipeline ##
# Fastqc
if [[ ! -e "${reads/.f*/_fastqc.zip}" ]] && [[ ! -z $fastQC ]]
then
  echo "Fastqc pretrimming"
  fastqc -t $threads "$reads" -o "$out_dir"
fi

# trim reads
qc_trim_SE $threads "$reads" "$trim_fasta" "$min_len" "$out_dir"

# align reads
bowtie --threads $threads --seed 1987 -x "$ref_rem" -q "$reads" -a -v $max_missmatch \
  --un "${reads}_cleaned.fq" > /dev/null
  # > "${nme}_rRNA.sam"
# samtools view -bS "${nme}_rRNA.sam" | samtools sort -@ $threads -O "bam" -T "working" -o "${nme}_rRNA.bam" 
# rm "${nme}_rRNA.sam"
# samtools index "${nme}_rRNA.bam"
# samtools flagstat "${nme}_rRNA.bam"

reads="${reads}_cleaned.fq"
bowtie --threads $threads -S --seed 1987 -x "$ref" -q "$reads" -a -v $max_missmatch "${out_dir}/${nme}.sam"
samtools view -bS "${out_dir}/${nme}.sam" | samtools sort -@ $threads -O "bam" -T "working" -o "${out_dir}/${nme}.bam"
samtools index "${out_dir}/${nme}.bam"
rm "${out_dir}/${nme}.sam"
samtools flagstat "${out_dir}/${nme}.bam"


# count alignments
echo "HtSeq started"
htseq-count --nprocesses $threads --type "gene" --idattr "Name" --order "name" --mode "union" --stranded "$strand" \
  --minaqual 10 --nonunique none -f bam "${out_dir}/${nme}.bam" "$gtf" --counts_output "${out_dir}/${nme}_HTSeq.gene.counts.tsv"

# featureCounts -F "GTF" -d 30 -s "$stran_fc" -t "gene" -g "Name" -O -Q 5 --ignoreDup -T $5 -a "$4" "$fCount" -o "${3/.bam/.featCount.counts}" "$3"

htseq-count --nprocesses $threads --type "CDS" --idattr "Name" --order "name" --mode "union" --stranded "$strand" \
  --minaqual 10 --nonunique none -f bam "${out_dir}/${nme}.bam" "$gtf" --counts_output "${out_dir}/${nme}_HTSeq.CDS.counts.tsv"



# 5′ mapped sites of RPFs
echo "psite started"
psite "${ref/.f*/_rois.txt}" "${nme}_riboprofile" \
  --min_length $min_len \
  --max_length $max_len \
  --require_upstream \
  --count_files "${out_dir}/${nme}.bam"


rm "${reads}"
##############        











