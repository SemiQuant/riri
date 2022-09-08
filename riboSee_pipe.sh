
TODO
table of metagenes to do, each colum a new list of genes
make the needed files, and report the number of nege found and not found in the original file






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
    -ms|--mask = mask stable RNAs in reference instead of prealigning?
    -u|--umi = UMI sequence if present #GNNNNNNNNGACTGGAGTTCAGACGTGTGCTCTTCCGA


    -p|--prime (defult = 3) plastid three or 5 prime
    -os|--offset plastid offset (defult = 14)
    -d|--downstream plastid downstream (defult = 100)
    -l|--landmark plastid landmark (defult = cds_start)
    -c|--codon_buffer plastid codon_buffer (defult = 5)
    -no|--normalize_over plastid normalize_over (defult = '30 200')
    -m|--min_counts plastid normalize_over (defult = 20)
    
    -pi|--plastid_input_extras A tsv file where each column is a list of genes of intrest, with the first entry the name of the list
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
  -ms|--mask)
  msk_or_rem="mask"
  ;;
  -u|--umi)
  umi="$2"
  ;;
  -p|--prime)
  plastid_prime="$2"
  ;;
  -os|--offset)
  offset="$2"
  ;;
  -d|--downstream)
  downstream="$2"
  ;;
  -l|--landmark)
  landmark="$2"
  ;;
 -c|--codon_buffer)
  codon_buffer="$2"
  ;;
 -no|--normalize_over)
  normalize_over="$2"
  ;;
 -m|--min_counts)
  min_counts="$2"
  ;;
  -pi|--plastid_input_extras)
  plastid_in="$2"
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
        -o "${out_fq}" "$2"
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
        fastqc -t $1 "${out_fq}" -o "${5}/"
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
  msk_or_rem="${msk_or_rem:-'Nomask'}"
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

  if [[ $plastid_prime == "5" ]]
  then
      plastid_prime = '--fiveprime '
  else
      plastid_prime = '--threeprime '
  fi

  offset="${offset:-14}"
  downstream="${downstream:-100}"
  landmark="${landmark:-cds_start}"

  codon_buffer="${codon_buffer:-5}"
  normalize_over="${normalize_over:-'30 200'}"
  min_counts="${min_counts:-20}"

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
    bowtie-build --quiet --threads $threads -f "$ref" "$ref"
  fi

  ref_rem="${ref/.f*/_rRNAs.fasta}"
  if [ ! -f "${ref_rem}.1.ebwt" ] 
  then
    # split rRNAs an tRNAs
    awk '{if ($3 == "tRNA" || $3 == "rRNA") print $0;}' "$gtf" > "${gtf/.g*/_rRNA.gtf}"
    if [[ $msk_or_rem == "mask" ]]
    then
      if [ ! -f "${ref/.f*/_rRNAsMasked.fasta}.1.ebwt" ]
      then
        bedtools maskfasta -fi "$ref" -bed "${gtf/.g*/_rRNA.gtf}" -fo "${ref/.f*/_rRNAsMasked.fasta}"
        bowtie-build --quiet --threads $threads -f "${ref/.f*/_rRNAsMasked.fasta}" "${ref/.f*/_rRNAsMasked.fasta}"
      fi
      ref="${ref/.f*/_rRNAsMasked.fasta}"
    else
      # I have no idea why things are getting thourgh, I've tried all the below, all have similar results
      bedtools getfasta -fi "$ref" -bed "${gtf/.g*/_rRNA.gtf}" > "$ref_rem"
      bowtie-build --quiet --threads $threads -f "$ref_rem" "$ref_rem"
      # bowtie cant handel multifasta, so just one per line
      # echo ">stableRNAs" > "$ref_rem"
      # bedtools getfasta -fi "$ref" -bed "${gtf/.g*/_rRNA.gtf}" | grep -v ">" | tr --delete '\n' >> "$ref_rem"
      # bowtie-build --quiet --threads $threads -f "$ref_rem" "$ref_rem"
      # stable_rnas="$(bedtools getfasta -fi "$ref" -bed "${gtf/.g*/_rRNA.gtf}" | grep -v ">" | tr '\n' ',')"
      # bowtie-build --quiet --threads $threads \
      #   -c "$stable_rnas" \
      #   "$ref_rem"
    fi
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
          export gtf_ps
        fi
      fi
          # run python script to make understandable gtf_tmp
          # this is only bcause bacterial annotations in file
          # python3 "${Script_dir}/gtf_primer.py" --gtf_in "$gtf_ps" > /dev/null 2>&1
    fi
        metagene generate "${ref/.f*/}" \
          --landmark $landmark \
          --downstream $downstream \
          --annotation_files "$gtf_ps"
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
  if [[ $msk_or_rem != "mask" ]]
  then
    bowtie --threads $threads --seed 1987 -x "$ref_rem" -q "$reads" -a -v $max_missmatch \
      --un "${reads}_cleaned.fq" > /dev/null 
    reads="${reads}_cleaned.fq"
  fi


  if [ ! -z "$umi" ]
  then
    # GNNNNNNNNGACTGGAGTTCAGACGTGTGCTCTTCCGA
    
    # umi_tools extract --stdin="$reads" \
    #   --bc-pattern=NNNNNNNNNCCCCCCCCCCCCCCCCCCCCCCCCCCCCC \
    #   --log="${nme}_processed.log" \
    #   --stdout "${nme}_UMI.fastq.gz"

      umi_tools extract --stdin="$reads" \
        --extract-method=regex \
        --bc-pattern='.*(?P<umi_1>.{9})(?P<discard_1>TCGGAAGAGCACACGTCTGAACTCCAGTC){s<=2}.*' \
        --log="${nme}_processed.log" \
        --stdout "${nme}_UMI.fastq.gz"

      reads="${nme}_UMI.fastq.gz"
  fi


  bowtie --threads $threads -S --seed 1987 -x "$ref" -q "$reads" -m 5 --best --strata -a -v $max_missmatch "${out_dir}/${nme}.sam"
  samtools view -bS "${out_dir}/${nme}.sam" | samtools sort -@ $threads -O "bam" -T "working" -o "${out_dir}/${nme}.bam"
  bam="${out_dir}/${nme}.bam"

  if [ ! -z "$umi" ]
  then
      umi_tools dedup -I "$bam" --output-stats="${out_dir}/${nme}.deduplicated" -S "${out_dir}/${nme}.deduplicated.bam"
      # sort?
      samtools sort -@ $threads -O "bam" -T "working" -o "${out_dir}/${nme}.deduplicated.sorted.bam" "${out_dir}/${nme}.deduplicated.bam"
      rm "${out_dir}/${nme}.deduplicated.bam"
      bam="${out_dir}/${nme}.deduplicated.sorted.bam"
  fi

  samtools index "$bam"
  rm "${out_dir}/${nme}.sam"
  samtools flagstat "$bam" > "${out_dir}/${nme}.flagstat"

  # count alignments
  echo "HtSeq started"
  htseq-count --nprocesses $threads --type "gene" --idattr "Name" --order "name" --mode "union" --stranded "$strand" \
    --minaqual 10 --nonunique none -f bam "$bam" "$gtf" --counts_output "${out_dir}/${nme}_HTSeq.gene.counts.tsv"


  htseq-count --nprocesses $threads --type "CDS" --idattr "Name" --order "name" --mode "union" --stranded "$strand" \
    --minaqual 10 --nonunique none -f bam "$bam" "$gtf" --counts_output "${out_dir}/${nme}_HTSeq.CDS.counts.tsv"



  if [[ $strand == "reverse" ]]
  then
      stran_fc=2
  elif [[ $strand == "yes" ]]
  then
      stran_fc=1
  else
      stran_fc=0
  fi

  featureCounts -F "GTF" -d 30 -s "$stran_fc" -t "gene" -g "Name" -O -Q 5 \
    --ignoreDup -T $threads -a "$gtf" "$fCount" -o "${out_dir}/${nme}_featCount.counts" "$bam"


  # 5′ mapped sites of RPFs
  # Not working
  # # generate metagene `roi` file
  # if [ ! -f "${ref/.f*/_rois.txt}" ]
  # then
  #     if [ ! -f "${gtf_ps}_5p.gtf" ]
  #     then
  #         if [[ "${gtf##*.}" == "gff" ]]
  #         then
  #           if [ -f "${gtf/.gff/.gtf}" ]
  #           then
  #             gtf_ps="${gtf/.gff/.gtf}"
  #           else
  #             gffread "$gtf" -T -o "${gtf/.gff/_tmp.gtf}"
  #             gtf_ps="${gtf/.gff/_tmp.gtf}"
  #             fi
  #           fi
  #         # run python script to make understandable gtf_tmp
  #         # this is only bcause bacterial annotations in file
  #         python3 "${Script_dir}/gtf_primer.py" --gtf_in "$gtf_ps" > /dev/null 2>&1
  #     fi
  #     metagene generate "${ref/.f*/}" \
  #       --landmark cds_start \
  #       --annotation_files "${gtf_ps}_5p.gtf"
  # fi

  # echo "psite started"
  # psite "${ref/.f*/_rois.txt}" "${nme}_riboprofile" \
  #   --min_length $min_len \
  #   --max_length $max_len \
  #   --require_upstream \
  #   --count_files "${out_dir}/${nme}.bam"


  # metagene generate NC_000962 \
  #     --landmark cds_start \
  #     --annotation_files NC_000962.gtf
  phase_by_size "${ref/.f*/_rois.txt}" "${nme}_phase_by_size" \
    --count_files "$bam" \
    "$plastid_prime" \
    --offset $offset \
    --codon_buffer "$codon_buffer" \
    --min_length $min_len \
    --max_length $max_len

  mkdir "${out_dir}/count_vectors"
  get_count_vectors --annotation_files "${ref/.f*/_rois.bed}" \
    --annotation_format BED \
    --count_files "$bam" \
    "$plastid_prime" \
    --offset $offset \
    --min_length $min_len \
    --max_length $max_len \
    "${out_dir}/count_vectors"

  zip -q -r "${out_dir}/count_vectors.zip" "${out_dir}/count_vectors"
  rm -r "${out_dir}/count_vectors"

  # '''
  # # from plastid import BAMGenomeArray, FivePrimeMapFactory, BED_Reader, Transcript
  # # transcripts = BED_Reader("/wynton/home/ernst/jdlim/riboseq/delete/NC_000962_rois.bed", return_type=Transcript)
  # # for transcript in transcripts:
  # #         print(transcripts)

  # import numpy
  # my_vector = numpy.loadtxt("/wynton/home/ernst/jdlim/riboseq/delete/count_vectors/Rv1669.txt")


  # # 30-codon sliding window average
  # window = numpy.ones(90).astype(float)/90.0
  # sliding_window_avg = numpy.convolve(my_vector,window,mode="valid")


  # # plot
  # import matplotlib.pyplot as plt

  # plt.plot(my_vector)
  # plt.plot(sliding_window_avg,label="30 codon average")
  # plt.xlabel("Position in transcript (5' to 3')")
  # plt.ylabel("Ribosome counts")

  # # add outlines at start & stop codons
  # # plt.axvline(my_transcript.cds_start, color="#999999",dashes=[3,2],zorder=-1)
  # # plt.axvline(my_transcript.cds_end, color="#999999",dashes=[3,2],zorder=-1)

  # plt.legend()
  # plt.show()
  # # plt.savefig( "test.png")
  # '''

  metagene count "${ref/.f*/_rois.txt}" "${nme}_count" \
    --count_files "$bam" \
    --offset $offset \
    "$plastid_prime" \
    --normalize_over "$normalize_over" \
    --min_counts "$min_counts" \
    --cmap Blues

  metagene chart "${nme}_counts_strt.png" \
    "${nme}_count_metagene_profile.txt" \
    --landmark "start codon"

  metagene chart "${nme}_counts_peak.png" \
    "${nme}_count_metagene_profile.txt" \
    --landmark "highest ribosome peak"

  if [[ $msk_or_rem != "mask" ]]
  then
    rm "$reads"
  fi

  ##############        


if [ ! -z "$gtf_ps" ]
then
    mkdir "${out_dir}/rois"
    cols=$(head -n1 $gtf_ps | grep -o "\t" | wc -l)
    for ColNum in {1..$cols}
    do
        nm="$(head -n1 $gtf_ps)"
        awk ‘{print $ColNum}’ "$gtf_ps" > "${out_dir}/rois/${nm}_plastid_roi.txt"
        grep -f "${nm}_plastid_roi.txt" > "${out_dir}/rois/${nm}_plastid_roi.gff"
        rm "${out_dir}/rois/${nm}_plastid_roi.txt"
        
        metagene generate "${out_dir}/rois/${nm}_plastid_roi" \
          --landmark $landmark \
          --downstream $downstream \
          --annotation_files "${out_dir}/rois/${nm}_plastid_roi.gff"
          
        metagene count "${out_dir}/rois/${nm}_plastid_roi_rois.txt}" "${nme}_count" \
          --count_files "$bam" \
          --offset $offset \
          "$plastid_prime" \
          --normalize_over "$normalize_over" \
          --min_counts "$min_counts" \
          --cmap Blues
      
        metagene chart "${out_dir}/rois/${nm}_plastid_roi_counts_strt.png" \
          "${out_dir}/rois/${nm}_plastid_roi_count_metagene_profile.txt" \
          --landmark "start codon"
      
        metagene chart "${out_dir}/rois/${nm}_plastid_roi_counts_peak.png" \
          "${out_dir}/rois/${nm}_plastid_roi_count_metagene_profile.txt" \
          --landmark "highest ribosome peak"
        
    done
fi

  # sort GTF file 
  # cat my_file.gtf | grep -v "#" | sort -k1,1 -k4,4n >my_file_sorted.gtf










