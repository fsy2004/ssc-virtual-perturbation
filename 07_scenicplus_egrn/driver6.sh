#!/bin/bash
# Driver v6 (post-biomart-403 fix). Annotation pre-downloaded+UCSC-fixed, download rule skipped.
# Runs snakemake (--cores 8) eGRN inference -> scplusmdata.h5mu -> extract. Disconnect-proof (setsid+nohup).
D=/data/ssc/scenicplus
DL=$D/logs/driver.log
source /data/ssc/miniconda3/etc/profile.d/conda.sh && conda activate scenicplus
cd $D
say(){ echo "[driver $(date '+%F %T')] $*" >> $DL; }
say "===DRIVER6 START (snakemake, annotation pre-downloaded UCSC-fixed, download rule skipped)==="
swapon --show 2>/dev/null | grep -q /data/swapfile || swapon /data/swapfile 2>/dev/null
free -g | awk 'NR==2{print "[driver] mem avail="$7}' >> $DL
SM=$D/scplus_pipeline/Snakemake
[ -s "$SM/genome_annotation.tsv" ] || { say "FAIL: genome_annotation.tsv missing"; exit 1; }
[ -s "$SM/chromsizes.tsv" ]        || { say "FAIL: chromsizes.tsv missing"; exit 1; }
say "genome_annotation.tsv + chromsizes.tsv present (UCSC style)"
cd $SM
snakemake --unlock --cores 1 >/dev/null 2>&1   # release any stale lock from the killed run
say "START snakemake (--cores 8)"
nice -n 10 snakemake --cores 8 --rerun-incomplete --keep-going >> $D/logs/stage_snakemake.log 2>&1
rc=$?; cd $D
MD=$(ls $SM/scplusmdata.h5mu 2>/dev/null | head -1)
if [ $rc -ne 0 ] || [ -z "$MD" ]; then say "FAIL snakemake (rc=$rc, mdata='$MD') — see logs/stage_snakemake.log"; exit 1; fi
say "DONE snakemake: $MD"
say "START extract"
nice -n 10 python -u run_scenicplus_pipeline.py extract >> $D/logs/stage_extract.log 2>&1
rc=$?
if [ $rc -ne 0 ] || ! grep -q "STAGE extract DONE" $D/logs/stage_extract.log; then
  say "FAIL extract (rc=$rc) — see logs/stage_extract.log"; exit 1; fi
say "DONE extract"
say "===DRIVER6 ALL DONE==="
