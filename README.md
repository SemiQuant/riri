# riri

## Pipelines to process bacterial PE RNAseq and SE RiboSeq
IInfo about pipeines will go here

I have written little to no error handling, so check logs etc.

## Download singularity container
./RNAseeker_pipe.sh --container


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
| -tm|--trim\_fasta | Path to adapter and linkers multi fasta | ${Script\_dir}/references/adapts.fasta |



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
