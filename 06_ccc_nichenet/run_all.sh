#!/bin/bash
# Disconnect-proof orchestrator for niche->myofibroblast Notch CCC (LIANA + NicheNet).
# Runs gently (nice) with capped threads so SCENIC+ keeps its cores.
OUT=/data/ssc/ccc_nichenet
LOG=$OUT/logs
mkdir -p "$LOG"
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
PY=/data/ssc/miniconda3/envs/npc/bin/python
RS=/data/ssc/miniconda3/envs/r44/bin/Rscript
ts(){ date '+%F %T'; }
echo "START $(ts)" > "$OUT/PIPELINE_STATUS"

run(){ # name cmd...
  local name="$1"; shift
  echo "[$(ts)] >>> $name" | tee -a "$OUT/PIPELINE_STATUS"
  nice -n 10 "$@" > "$LOG/${name}.log" 2>&1
  local rc=$?
  echo "[$(ts)] <<< $name rc=$rc" | tee -a "$OUT/PIPELINE_STATUS"
  if [ $rc -ne 0 ]; then echo "FAILED $name rc=$rc $(ts)" >> "$OUT/PIPELINE_STATUS"; exit $rc; fi
}

run 00_atlas_and_liana $PY $OUT/00_atlas_and_liana.py
run 02_nichenet        $RS $OUT/02_nichenet.R
run 03_report_figs     $PY $OUT/03_report_figs.py
echo "DONE $(ts)" >> "$OUT/PIPELINE_STATUS"
