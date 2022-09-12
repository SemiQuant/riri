# riri

## Pipelines to process bacterial PE RNAseq and SE RiboSeq
RiboSee is used to align singe end Ribosomal profiling reads (RiboSeq) to a reference, and output alignments, counts, psite determination, etc.
RNAseeker is used to align pair-end sequencing to a reference, and output alignments, counts and statistics.

Secondary scripts can be used for analyses.
- Differential RNAseq
- Psite statistical anlayses
- Etc

I have written little to no error handling, so check logs etc.
There is also something strange happenign with the prealignemnts to the stable RNAs and many are getting through, so I masked these in the reference just incase.

If files dont execute then do this
`chmod +x RNAseeker_pipe.sh riboSee_pipe.sh gtf_primer.py singularity_continer_setup.sh`


# TODO
-Fix error with qualimap? (Error while calculating counts! Failed to detect annotations file format.


## Download scripts
`git clone --recursive https://github.com/SemiQuant/riri.git`

## Download singularity container
`./RNAseeker_pipe.sh --container`


## RiboSee

| Flag | Description | Defaults |
| --- | --- | --- |
| -t|--threads | Number of threads to use | NA |
| -g|--genome\_reference | Full path to reference genome | ${Script\_dir}/references/NC_000962.fasta |
| -gtf|--GTF\_reference | Full path to reference annotations | ${Script\_dir}/references/NC_000962.gff |
| -rd|--read\_dir | Full path to loaction of read file | NA |
| -r|--reads | Read name | NA |
| -o|--out\_dir | Full path to output directory | NA |
| -n|--name | Sample name | NA |
| -s|--strand | Stranded sequecning (yes|no|reverse) | no |
| -sd|--script\_directory | Path to script location | $(dirname &quot;$0&quot;) |
| -dl|--container | Downloaad singularity container and exit | NA |
| -mt|--get\_metrics | Compile report of metric after run, path to folder of results | NA |
| -fq|--fastQC | Perform fastQC anlsysis | NA |
| -mm|--max\_missmatch | Maximum missmatches for alignment | 2 |
| -mn|--min\_len | Mininum read length | 24 |
| -mx|--max\_len | Maximum read length | 36 |
| -tm|--trim\_fasta | Path to adapter and linkers multi fasta, uses Trimmomatic | ${Script\_dir}/references/adapts.fasta |
| -ca|--cut\adapt | Adapter sequence to cut (e.g., CTGTAGGCACCATCAAT); Overwrites trim_fasta and uses CutAdapter | NA |
| -ms|--mask | mask stable RNAs in reference instead of prealigning to them? | NA |
| -u|--umi | UMI sequence if present #GNNNNNNNNGACTGGAGTTCAGACGTGTGCTCTTCCGA | NA |
| -p|--prime | (defult = 3) plastid three or 5 prime | NA |
| -os|--offset | plastid offset (defult = 14) | NA |
| -d|--downstream | plastid downstream (defult = 100) | NA |
| -l|--landmark | plastid landmark (defult = cds_start) | NA |
| -c|--codon_buffer | plastid codon_buffer (defult = 5) | NA |
| -no|--normalize_over | plastid normalize_over (defult = '30 200') | NA |
| -m|--min_counts | plastid normalize_over (defult = 20) | NA |
| -pi|--plastid_input_extras | A tsv file where each column is a list of genes of intrest, with the first entry the name of the list | NA |


### Example run
also see "wynton_slurm_wrapper.sge"

```
out_dir="/wynton/home/ribSeq"
container="/wynton/home/riri_v0.1.sif"
script_dir="/wynton/home/riri"
read_dir="/wynton/home/fastq"
nm="file_name"

mkdir -p "$out_dir"
cd "$out_dir"

singularity exec "$container" \
  "${script_dir}/riboSee_pipe.sh" \
  --threads 8 \
  --genome_reference "${script_dir}/references/NC_000962_rRNAsMasked.fasta" \
  --GTF_reference "${script_dir}/references/NC_000962.gff" \
  --reads "${read_dir}/${nm}_L2_1.fq.gz" \
  --out_dir "$out_dir" \
  --name "$nm" \
  --strand "reverse" \
  --script_directory "${script_dir}" \
  --fastQC \
  --max_missmatch 2 \
  --min_len 24 \
  --max_len 36 \
  --trim_fasta "${script_dir}/references/adapts.fasta"
```


## RNAseeker

| Flag | Description | Defaults |
| --- | --- | --- |
| -r|--ref | Full path to reference genome | ${Script\_dir}/references/NC_000962.fasta |
| -t|--threads | Number of threads to use |
| -g|--gtf | Full path to reference annotations | ${Script\_dir}/references/NC_000962.gff |
| -r1|--read1 | Full path to loaction of read1 file | NA |
| -r2|--read2 | Full path to loaction of read2 file | NA |
| -n|--name | Sample name | NA |
| -o|--out\_dir | Full path to output directory | NA |
| -m|--ram | Ram | 2\*threads |
| -s|--strand | Stranded sequecning (yes|no|reverse) | no |
| -rR|--remove\_rRNA | Remove rRNA from annotation file | NA |
| -sd|--script\_directory | Path to script location | $(dirname &quot;$0&quot;) |
| -dl|--container | Downloaad singularity container and exit | NA |
| -mt|--get\_metrics | Compile report of metric after run, path to folder of results | NA |
| -tr|--trim\_metrics | trim reads? | NA |
| -a|--adapters | Path to adapter and linkers multi fasta | ${Script\_dir}/references/adapts.fasta |
| -fq|--fastQC | Perform fastQC anlsysis | NA |


### Example run
also see "wynton_slurm_wrapper.sge"

```
out_dir="/wynton/home/RNAseq"
container="/wynton/home/riri_v0.1.sif"
script_dir="/wynton/home/riri"
read_dir="/wynton/home/fastq"
nm="file_name"

mkdir -p "$out_dir"
cd "$out_dir"

singularity exec "$container" \
  "${script_dir}/RNAseeker_pipe.sh" \
    --ref "${script_dir}/references/NC_000962.fasta" \
    --threads 8 \
    --gtf "${script_dir}/references/NC_000962.gff" \
    --read1 "${read_dir}/${nm}_L2_1.fq.gz" \
    --read2 "${read_dir}/${nm}_L2_2.fq.gz" \
    --name ${nm} \
    --out_dir "$out_dir" \
    --adapters "${script_dir}/references/adapts.fasta" \
    --strand "reverse" \
    --trim \
    --remove_rRNA \
    --fastQC \
    --keep_unpaired \
    --script_directory "${script_dir}"
```
