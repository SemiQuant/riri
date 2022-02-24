#!/bin/bash
usage () { #echo -e to read \n
  echo "
Usage Options
  Mode 1:
    -dl|--container = download the singularity container to this path and exit
  
  Mode 2:
    -mt|--get_metrics = supply a dir and get metrics for all analyses in that dir, all other options will be ignored if this is non-empyt
  
  Mode 3: 
    -r|--ref = genome fasta
    -t|--threads
    -g|--gtf = gtf file
    -r1|--read1
    -r2|--read2 = the second fastq file if PE reads
    -n|--name = output name
    -o|--out_dir
    -m|--ram
    -a|--adapters multifasta of adapters to clip
    -s|--strand stranded library (yes|no|reverse)
    -tr|--trim = trim reads?
    # -tm|--trim_min (if trim flag selected)
    -k|--keep_unpaired (if trim flag selected)
    -rR|--remove_rRNA = remove rRNA from annotation file
    -sd|--script_directory
    -fq|--fastQC = run fastqc?
  "
}

if [ $# == 0 ]
then
    usage
    exit 1
fi


declare_globals () {
    # if same thing twice will take second one
    while [[ "$#" -gt 0 ]]
    do
        case $1 in
        -r|--ref)
        ref="$2"
        ;;
        -t|--threads)
        threads="$2"
        ;;
        -g|--gtf)
        gt="$2"
        ;;
        -r1|--read1)
        read1="$2"
        ;;
        -r2|--read2)
        read2="$2"
        ;;
        -n|--name)
        name="$2"
        ;;
        -o|--out_dir)
        out_dir="$2"
        ;;
        -m|--ram)
        ram="$2"
        ;;
        -a|--adapters)
        adapters="$2"
        ;;
        -s|--strand)
        strand="$2"
        ;;
        # -tm|--trim_min)
        # trim_min="$2"
        # ;;
        -rR|--remove_rRNA)
        rRNA="Y"
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
        -tr|--trim)
        trim="Y"
        ;;
    esac
        shift
    done
}


extract () {
    if [ -f $1 ] ; then
      case $1 in
        *.tar.bz2)   tar xjf $1     ;;
        *.tar.gz)    tar xzf $1     ;;
        *.bz2)       bunzip2 $1     ;;
        *.rar)       unrar e $1     ;;
        *.gz)        gunzip $1      ;;
        *.tar)       tar xf $1      ;;
        *.tbz2)      tar xjf $1     ;;
        *.tgz)       tar xzf $1     ;;
        *.zip)       unzip $1       ;;
        *.Z)         uncompress $1  ;;
        *.7z)        7z x $1        ;;
        *)  echo "'$1' is not zipped" 2>/dev/null ;;
         esac
    fi
}


qc_trim_PE () {
    #FastQC pre
    if [[ ! -e "${3}/${1/.f*/_fastqc.zip}" ]] || [[ ! -e "${3}/${2/.f*/_fastqc.zip}" ]]
    then
        if [[ ! -z $fastQC ]]
        then
            fastqc -t $5 "$1" -o "$3"
            fastqc -t $5 "$2" -o "$3"
        fi
    fi
    
    #Trim Reads
    echo "trimming started $1 $2"
    if [[ ! -z $trim ]]
    then
        # if [[ -e "${1/f*/forward.fq.gz}" ]]
        read1="${3}/$(basename ${1/.f*/_forward.fq.gz})"
        read2="${3}/$(basename ${2/.f*/_reverse.fq.gz})"
        if [[ -e "$read1" ]] && [[ -e "$read2" ]]
        then
            echo "Found ${1/.f*/}"
            export read1
            export read2
        else
            java -jar "$TRIM" PE -phred33 \
              -threads $5 \
              "$1" "$2" \
              "${1/.f*/_forward_paired.fq.gz}" "${1/.f*/_forward_unpaired.fq.gz}" \
              "${2/.f*/_reverse_paired.fq.gz}" "${2/.f*/_reverse_unpaired.fq.gz}" \
              ILLUMINACLIP:"$6":2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:$7
        
            # mv "${1/.f*/_forward.fq.gz}" "${2/.f*/_reverse.fq.gz}" "$3"
            mv "${1/.f*/_forward_paired.fq.gz}" "$read1"
            mv "${2/.f*/_reverse_paired.fq.gz}" "$read2"
            
            if [[ ! -z $keep_unpaired ]]
            then
                #as we also want unpaired reads so..
                # woudl be quicker if you do a cp and then a cat
                # cat "${1/.f*/_forward_paired.fq.gz}" "${1/.f*/_forward_unpaired.fq.gz}" > "$read1"
                # cat "${2/.f*/_reverse_paired.fq.gz}" "${2/.f*/_reverse_unpaired.fq.gz}" > "$read2"
                cat "${1/.f*/_forward_unpaired.fq.gz}" >> "$read1"
                cat "${2/.f*/_reverse_unpaired.fq.gz}" >> "$read2"
                
                # just to make sure as sometimes it lets one through if merging the paired and unpaired
                cutadapt --cores=$5 --quality-cutoff 10,10 --minimum-length $7 -o "${read1}.tmp.gz" "$read1"
                mv "${read1}.tmp.gz" "$read1"
                cutadapt --cores=$5  --quality-cutoff 10,10 --minimum-length $7 -o "${read2}.tmp.gz" "$read2"
                mv "${read2}.tmp.gz" "$read2"
            fi
            
            mv "${1/.f*/_forward_unpaired.fq.gz}" "${2/.f*/_reverse_unpaired.fq.gz}" "${3}"
            
            #FastQC post
            if [[ ! -z $fastQC ]]
            then
                fastqc -t $5 "$read1" -o "$3"
                fastqc -t $5 "$read2" -o "$3"
            fi
        fi
        export read1
        export read2
    fi
    echo "trimming completed"
}



BOWTIE_index () {
  if [ ! -e "${1/.f*/}.1.bt2" ] #check if indexed alread #${1/.f*/.1.bt2}
  then
      if [[ $1 == *.gz ]]
      then
          gunzip $1
          bowtie2-build --threads $2 ${1/.gz/} ${1/.f*/} #${1/.f*/} #$(printf $1 | cut -f 1 -d '.')
      else
          bowtie2-build --threads $2 $1 ${1/.f*/} #${1/.f*/} #$(printf $1 | cut -f 1 -d '.')
      fi
  else
      echo "Found bowtie2 index for $1?"
  fi
}


BOWTIE_alignerPE () {
    echo "BOWTIE alignment started $3"
    
    # BOWTIE_alignerPE "$read1_unaligned" "$threads" "$g2" "$out_dir" "$name" "$ram" "$read2_unaligned"
    local ref="${3/.f*/}"
    local gen=$(basename $ref)
    local out_f="${4}/${5}.$(printf $gen | cut -f 1 -d '.').sam"
    local out_f_bam="${out_f/sam/bam}"
    local R1="$1"
    local R2="$7"
    local thread=$2
    
    if [[ -e "$out_f_bam" ]]
    then
        echo "Found $out_f_bam"
        export read1_unaligned="${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz" 
        export read2_unaligned="${4}/${5}_${gen}_Unmapped.out.mate2.fastq.gz"
    else
        # bowtie2 --n-ceil L,0,0.05 --score-min L,1,-0.6 -p "$2" -x ${3/.f*/}  -1 "$1" -2 "$7" -S "$out_f" --un-gz ${4} --un-conc-gz ${4}
         bowtie2 \
            --dovetail \
            --local \
            --minins 0 \
            -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 \
            -p "$thread" -x "$ref" -1 "$R1" -2 "$R2" -S "$out_f" --un-gz "${4}" --un-conc-gz "${4}"

    export read1_unaligned="${4}/${5}_${gen}_Unmapped.out.mate1.fastq.gz" 
    export read2_unaligned="${4}/${5}_${gen}_Unmapped.out.mate2.fastq.gz"

    mv "un-conc-mate.1" "$read1_unaligned" 
    mv "un-conc-mate.2" "$read2_unaligned"
    
    java -jar "$PICARD" SortSam \
      I="$out_f" \
      O="$out_f_bam"\
      SORT_ORDER=queryname \
      VALIDATION_STRINGENCY=LENIENT

    rm "$out_f"
    
    fi
    echo "BOWTIE alignment completed"
}




do_calcs () {
    if [[ $strand == "reverse" ]]
    then
        stran_fc=2
        stran_qm="strand-specific-reverse"
        # LT=
    elif [[ $strand == "yes" ]]
    then
        stran_fc=1
        stran_qm="strand-specific-forward"
    else
        stran_fc=0
        stran_qm="non-strand-specific"
    fi
    
    
    fCount='-p' #this sets it to PE
    
    samtools flagstat "$3" > "${3/.bam/_flagstat.txt}"



    echo "Started htseq-count $(basename $2)"
        htseq-count --type "gene" --idattr "Name" --order "name" --mode "union" --stranded "$strand" -a 5 --nonunique all -f bam "$3" "$4" > "${3/bam/HTSeq.counts}" #
        
        # or gene? - let user input type to count
        if [[ ! -z $feat ]]
        then
            echo "Started featureCounts $(basename $2)"
            featureCounts -F "GTF" -d 30 -s "$stran_fc" -t "gene" -g "Name" -O -Q 5 --ignoreDup -T $5 -a "$4" "$fCount" -o "${3/.bam/.featCount.counts}" "$3"
        fi

    echo "Counts completed"
    
    
    #can also do qualimap
    # export PATH=/users/bi/jlimberis/bin/qualimap_v2.2.1:$PATH
    if [[ ! -z $qualimap ]]
    then
        #qualimap rnaseq only works with gtf file??
        gtf="$4"
        if [[ "${gtf##*.}" == "gff" ]]
        then
            gffread "$gtf" -T -o "${gtf/.gff/_tmp.gtf}"
            local gtf="${gtf/.gff/_tmp.gtf}"
        fi
            # samtools sort -n -@ $5 -o "${3/bam/coord.bam}" "$3"
            # --sorted is giving issues.. have to let it do it? this take much more time
            echo "Started qualimap rnaseq $(basename $2)"
            qualimap rnaseq --paired -p "$stran_qm" -bam "$3" -gtf "$gtf" -outdir "${3/.bam/_qualimap}" 1>/dev/null
            # --sorted
            
            echo "Started qualimap comp-counts $(basename $2)"
            sed 's/exon/CDS/g' "$gtf" > "${gtf}.tmp"
            qualimap comp-counts -bam "$3" -gtf "${gtf}.tmp" -id "gene_name" -type "CDS" -s -out "${3/.bam/_counts.html}" -p "$stran_qm" --paired 1>/dev/null
            # --sorted
            
            rm "${gtf}.tmp"

    fi
}


Multi_met_pic () {
    # $1 = gff
    # $2 = threads
    # $3 bam
    # $4 name
    # $5 ref
    # $6 strand
    
    #PICARD requires coordinate sorted bam file
    if [[ ! -e "${3/bam/coord.bam}" ]]
    then
        samtools sort -@ $2 -o "${3/bam/coord.bam}" "$3"
    fi
    
    java -jar "$PICARD" CollectMultipleMetrics \
        I="${3/bam/coord.bam}" \
        O="${4}_multiple_metrics" \
        R="$5"
      

    local gen=$(basename $1)
    local gen=${gen/.g*/}
    local gtf_1"=$1"
    
    
    if [[ ! -e "${gtf_1%.*}.rRNA.gtf" ]]
    then
        cat $gtf_1 | grep "gbkey=rRNA" | grep -v "ribosomal RNA protein" > "${gtf_1%.*}.rRNA.gff" # All rRNAs
        if [[ "${gtf_1##*.}" == "gff" ]]
        then
            gffread "${gtf_1%.*}.rRNA.gff" -T -o "${gtf_1%.*}.rRNA.gtf"
        fi
        cat "${gtf_1%.*}.rRNA.gtf" | sed 's/\texon\t/\tgene\t/g' > tmp_rRNA # Add gene lines
        cat "${gtf_1%.*}.rRNA.gtf" | sed 's/\texon\t/\ttranscript\t/g' >> tmp_rRNA # Add transcript lines
        cat tmp_rRNA >> "${gtf_1%.*}.rRNA.gtf"
        rm tmp_rRNA
        sed -i -e 's/$/ gene_biotype \"rRNA\"; transcript_biotype \"rRNA\"; gene_source \"ncbi\"; transcript_source \"ncbi\";/' "${gtf_1%.*}.rRNA.gtf" # Add to end of each line
    fi
    
    if [[ ! -e "${gtf_1%.*}.refFlat.txt" ]]
    then
        gtfToGenePred -genePredExt "${gtf_1%.*}.rRNA.gtf" refFlat.txt.tmp
        paste <(cut -f 12 refFlat.txt.tmp) <(cut -f 1-10 refFlat.txt.tmp) > "${gtf_1%.*}.refFlat.txt" # We have to add the gene name to the refFlat https://www.biostars.org/p/120145/; #cat ${GTF%.*}.refFlat.txt | awk 'BEGIN{FS="\t";} {OFS="\t";} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF%.*}.refFlat.txt # Same as above but written in a different way
        rm refFlat.txt.tmp
        
        gffread "$gtf_1" -T -o "${gtf_1%.*}.all.rRNA.gtf"
        gtfToGenePred -genePredExt "${gtf_1%.*}.all.rRNA.gtf" refFlat2.txt.tmp
        paste <(cut -f 12 refFlat2.txt.tmp) <(cut -f 1-10 refFlat2.txt.tmp) > "${gtf_1%.*}.all.refFlat.txt" # We have to add the gene name to the refFlat https://www.biostars.org/p/120145/; #cat ${GTF%.*}.refFlat.txt | awk 'BEGIN{FS="\t";} {OFS="\t";} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > ${GTF%.*}.refFlat.txt # Same as above but written in a different way
        rm refFlat2.txt.tmp
    fi
    
    
    # Intervals for rRNA transcripts.
    samtools view -H "$bam" > "${4}_rRNA.intervalListBody.txt"

    cat ${gtf_1%.*}.rRNA.gtf | awk '$3 == "transcript"' | \
      cut -f1,4,5,7,9 | \
      perl -lane '
          /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
          print join "\t", (@F[0,1,2,3], $1)
      ' | \
      sort -k1V -k2n -k3n \
      >> "${4}_rRNA.intervalListBody.txt"
      
    
    if [[ "$6" == "reverse" ]]
    then
        local strandP="SECOND_READ_TRANSCRIPTION_STRAND"
    elif [[ $strand == "yes" ]]
    then
        local strandP="FIRST_READ_TRANSCRIPTION_STRAND"
    else
        local strandP="NONE"
    fi
    
        
    java -jar "$PICARD" CollectRnaSeqMetrics \
        I="${3/bam/coord.bam}" \
        O="${4}_RNA_Metrics" \
        REF_FLAT="${gtf_1%.*}.all.refFlat.txt" \
        STRAND_SPECIFICITY="$strandP" \
        RIBOSOMAL_INTERVALS="${4}_rRNA.intervalListBody.txt" \
        CHART_OUTPUT="${4}.pdf"
        #ASSUME_SORTED=FALSE

    rm "${3/bam/coord.bam}"
}






############
# pipeline #
# setup variables
declare_globals "$@"

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


# Script_dir_tmp=$(dirname "$0")
Script_dir_tmp="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
Script_dir="${Script_dir:-$Script_dir_tmp}"
ram_def=$(expr $threads \* 2)
ram="${ram_in:-$ram_def}"
jav_ram=$(echo "scale=2; $ram*0.8" | bc)
export _JAVA_OPTIONS=-Xmx"${jav_ram%.*}G"
strand="${strand:-no}"
trim_min=16

trim_tmp="${Script_dir}/references/adapts.fasta"
trim_fasta="${trim_fasta:-trim_tmp}"
ref_tmp="${Script_dir}/references/Mycobacterium_tuberculosis_H37Rv_genome_v4.fasta"
ref="${ref:-ref_tmp}"
gtf_tmp="${Script_dir}/references/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff"
gtf="${gtf:-gtf_tmp}"


# PATHs in singularity container
TRIM=/usr/bin/Trimmomatic-0.39/trimmomatic-0.39.jar
adapters="${Script_dir}/references/adapts.fasta"
PICARD=/usr/bin/picard.jar

out_dir="${out_dir:-read_dir}"
mkdir "${out_dir}"
name="${name:-${read1/.f*/}}"

read1="$read_dir/$read1"
read2="$read_dir/$read2"

mkdir "${out_dir}/${name}"
out_dir="${out_dir}/${name}"
back_dir=${PWD}
cd "$out_dir"

#check if programs installed
command -v bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools but it's not installed. Aborting."; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo >&2 "I require bcftools but it's not installed. Aborting."; exit 1; }
command -v fastqc >/dev/null 2>&1 || { echo >&2 "I require fastqc but it's not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed. Aborting."; exit 1; }
command -v htseq-count >/dev/null 2>&1 || { echo >&2 "I require htseq but it's not installed. Aborting."; exit 1; }
command -v python >/dev/null 2>&1 || { echo >&2 "I require python2 but it's not installed. Aborting."; exit 1; }
command -v featureCounts >/dev/null 2>&1 || { echo >&2 "I require featureCounts but it's not installed. Aborting."; exit 1; }
command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "I require bowtie2 but it's not installed. Aborting."; exit 1; }
if [ ! -z $cut_adapt  ]; then command -v cutadapt >/dev/null 2>&1 || { echo >&2 "I require cutadapt but it's not installed. Aborting."; exit 1; }; fi

command -v python3 >/dev/null 2>&1 || { echo >&2 "I may require python3 but it's not installed."; }

if [ ! -f "$TRIM" ]; then echo "$TRIM not found!"; exit 1; fi
if [ ! -f "$PICARD" ]; then echo "$PICARD not found!"; exit 1; fi


# create index
BOWTIE_index "$ref" "$threads" "$gt"

# remove rRNA
if [[ ! -v $rRNA ]]
then
  # sed '/rRNA/d;/ribosomal RNA/d;/ribosomal/d' "$gt" > "$gt_no_rRNA"
  awk '{if ($3 != "rRNA") print $0;}' "$gt" > "${gt}_no_rRNA"
  export gt="${gt}_no_rRNA"
fi


# trim
if [[ -f $adapters ]]
then
  qc_trim_PE "$read1" "$read2" "$out_dir" $ram $threads "$adapters" $trim_min
fi


# alignments
read_length=$(zcat $read1 | head -n 10000 | awk '{if(NR%4==2) print length($1)}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')

BOWTIE_alignerPE "$read1" $threads "$ref" "$out_dir" "$name" $ram "$read2"
bam_file="${out_dir}/${name}.$(printf $(basename $ref) | cut -f 1 -d '.').bam"
do_calcs $out_dir $ref $bam_file $gt $threads $strand $read_length
Multi_met_pic "$gt" $threads "$bam_file" "$name" "$ref" "$strand"

# Cleanup dirs



cd "$back_dir"


