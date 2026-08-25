#!/bin/bash
GD="$HOME/.local/share/genomes/hg38"
mkdir -p "$GD"; cd "$GD"
if [ -s hg38.fa ] && [ -s hg38.fa.fai ]; then echo "GENOME_ALREADY_OK"; cat hg38.fa.sizes; exit 0; fi
URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
curl -fL --retry 12 --retry-delay 10 --max-time 4000 "$URL" \
 | zcat | awk '
   /^>/ { split($0,a," "); id=substr(a[1],2); keep=0; name="";
          if (id ~ /^([1-9]|1[0-9]|2[0-2])$/) {name="chr" id; keep=1}
          else if (id=="X"||id=="Y") {name="chr" id; keep=1}
          else if (id=="MT") {name="chrM"; keep=1}
          if (keep) print ">" name; next }
   { if (keep) print }' > hg38.fa
echo "awk_exit ${PIPESTATUS[*]}"
/data/ssc/miniconda3/envs/co2/bin/python -c "import pysam; pysam.faidx('$GD/hg38.fa')"
awk '{print $1"\t"$2}' hg38.fa.fai > hg38.fa.sizes
printf "hg38 = Ensembl GRCh38 release-110 primary_assembly, main chroms renamed to UCSC chr\n" > README.txt
echo "GENOME_SETUP_DONE"; echo "n_seqs:"; wc -l hg38.fa.fai; cat hg38.fa.sizes
