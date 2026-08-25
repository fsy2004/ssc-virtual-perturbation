#!/bin/bash
# Autonomous chain: waits STEP1(Cicero) -> STEP2(base-GRN parquet) -> gate(b) guard
#                   -> edit CFG(skin_base_grn+n_jobs) -> nohup KO re-run.
# Runs under nohup so it survives operator disconnects. Enforces the gate-(b)
# guard: KO launches ONLY if HES1 AND SMAD3 are present TFs with non-trivial targets.
cd /data/ssc/basegrn
CO2=/data/ssc/miniconda3/envs/co2/bin/python
ENGINE=/data/ssc/server_powered_2026/03_celloracle_engine.py
L=/data/ssc/basegrn/logs/chain.log
say(){ echo "[chain $(date '+%H:%M:%S')] $*" >> "$L"; }
: > "$L"; say "chain started"

# ---- 1. wait for STEP1 (Cicero) ----
STEP1_OK=0
for i in $(seq 1 170); do          # up to ~2.8h
  if grep -q "DONE_R" logs/build_r.log 2>/dev/null; then STEP1_OK=1; break; fi
  if grep -q "Execution halted" logs/build_r.log 2>/dev/null; then STEP1_OK=0; break; fi
  sleep 60
done
if [ "$STEP1_OK" != "1" ]; then say "STEP1 did NOT finish cleanly (no DONE_R). ABORT."; tail -8 logs/build_r.log >> "$L"; exit 1; fi
say "STEP1 done: $(grep DONE_R logs/build_r.log | tail -1)"

# ---- 2. STEP2 base-GRN build (motif scan) ----
say "STEP2 base-GRN build starting"
nice -n 5 "$CO2" 11_build_base_grn.py > logs/step2_basegrn.log 2>&1
say "STEP2 exit=$?"
tail -12 logs/step2_basegrn.log >> "$L"

# ---- 3. gate (b): HES1 AND SMAD3 present as TFs, parquet saved ----
SAVED=$(grep -c "^SAVED" logs/step2_basegrn.log)
HES1_OK=$(grep -c "TF HES1: present=True" logs/step2_basegrn.log)
SMAD3_OK=$(grep -c "TF SMAD3: present=True" logs/step2_basegrn.log)
if [ ! -s /data/ssc/basegrn/skin_base_GRN_dataframe.parquet ] || [ "$SAVED" -lt 1 ] || [ "$HES1_OK" -lt 1 ] || [ "$SMAD3_OK" -lt 1 ]; then
  say "GATE-B FAILED (parquet_saved=$SAVED HES1=$HES1_OK SMAD3=$SMAD3_OK). NOT launching KO. Needs debug."
  exit 2
fi
say "GATE-B PASSED: parquet built, HES1 & SMAD3 are source TFs. Proceeding to KO."

# ---- 4. edit CFG (ONLY skin_base_grn + n_jobs; scoring untouched) ----
cp -n "$ENGINE" "${ENGINE}.bak_skinarm"
grep -q 'skin_base_GRN_dataframe.parquet"' "$ENGINE" || \
  sed -i 's|"skin_base_grn": None,|"skin_base_grn": "/data/ssc/basegrn/skin_base_GRN_dataframe.parquet",|' "$ENGINE"
sed -i 's|"n_jobs": 18,|"n_jobs": 10,|' "$ENGINE"
say "CFG now: $(grep -n 'skin_base_grn\|\"n_jobs\"' "$ENGINE" | head -3 | tr '\n' ' ')"

# ---- 5. launch KO re-run (skin arm) ----
cd /data/ssc/server_powered_2026
nohup "$CO2" 03_celloracle_engine.py > /data/ssc/basegrn/logs/ko_skinarm.log 2>&1 &
say "KO re-run launched PID $! -> logs/ko_skinarm.log"
echo "CHAIN_LAUNCHED_KO" >> "$L"
