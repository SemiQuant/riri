# ririUMI
read=""
nme=""
umi_tools extract --stdin="$read" \
  --bc-pattern=GNNNNNNNNGACTGGAGTTCAGACGTGTGCTCTTCCGA \
  --log="${nme}_processed.log" \
  --stdout "${nme}_UMI.fastq.gz"



MAP


umi_tools dedup -I "${nme}.bam" --output-stats=deduplicated -S "${nme}.deduplicated.bam"

