#!/usr/bin/env python3
"""Independently validate a Tier-A low-energy +/-45 extension gate."""
from __future__ import annotations
import argparse,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument("--authorization",required=True);p.add_argument("--expected-rotor-id",required=True);p.add_argument("--output-dir",required=True);a=p.parse_args();ap=Path(a.authorization).resolve();x=json.loads(ap.read_text())
 if x["status"]!="pass_tier_a_plus_minus45_extension_authorized" or x["rotor_id"]!=a.expected_rotor_id or x["authorized_signed_step_indices"]!=[3,-3] or x["farther_extension_authorized"] is not False:raise ValueError("authorization structure mismatch")
 if not all(float(v)<=10.0 for v in x["terminal_relative_energies_kcal_mol"].values()):raise ValueError("low-energy trigger mismatch")
 for item in x["inputs"].values():
  if sha(item["path"])!=item["sha256"]:raise ValueError("input hash mismatch")
 out=Path(a.output_dir).resolve();out.mkdir(parents=True,exist_ok=False);r={"schema_version":"1.0","report_type":"o6u_tier_a_low_energy_extension_independent_validation","status":"pass_independent_tier_a_low_energy_extension","created_at_utc":datetime.now(timezone.utc).isoformat(),"authorization":{"path":str(ap),"sha256":sha(ap)},"rotor_id":a.expected_rotor_id,"validated_signed_step_indices":[3,-3],"farther_extension_authorized":False,"production_approved":False};f=out/"O6U_TIER_A_LOW_ENERGY_EXTENSION_INDEPENDENT_VALIDATION.json";f.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n");print(json.dumps({"status":r["status"],"report":str(f),"sha256":sha(f)},sort_keys=True))
if __name__=="__main__":main()
