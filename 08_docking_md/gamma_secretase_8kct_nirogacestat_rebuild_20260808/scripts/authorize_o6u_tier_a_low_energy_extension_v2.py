#!/usr/bin/env python3
"""Authorize exactly one Tier-A +/-45 extension from validated low-energy +/-30 points."""
from __future__ import annotations
import argparse,hashlib,json,math
from datetime import datetime,timezone
from pathlib import Path
K=627.5094740631
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(p):p=Path(p).resolve();return {"path":str(p),"size_bytes":p.stat().st_size,"sha256":sha(p)}
def main():
 p=argparse.ArgumentParser();p.add_argument("--adaptive-scope",required=True);p.add_argument("--positive-report",required=True);p.add_argument("--positive-validation",required=True);p.add_argument("--negative-report",required=True);p.add_argument("--negative-validation",required=True);p.add_argument("--rotor-id",required=True);p.add_argument("--frame-energy-hartree",type=float,required=True);p.add_argument("--output-dir",required=True);a=p.parse_args()
 paths={k:Path(v).resolve() for k,v in vars(a).items() if k.endswith("report") or k.endswith("validation") or k=="adaptive_scope"}
 scope=json.loads(paths["adaptive_scope"].read_text());pr=json.loads(paths["positive_report"].read_text());nr=json.loads(paths["negative_report"].read_text());pv=json.loads(paths["positive_validation"].read_text());nv=json.loads(paths["negative_validation"].read_text())
 if sha(paths["adaptive_scope"])!="519d587648c6458e23dca0e7cc9041a3a711feb6d793b55bd3efcdc532aaa49e" or scope["status"]!="pass_compound_specific_adaptive_scan_scope_authorized":raise ValueError("scope mismatch")
 for r,v,s in ((pr,pv,2),(nr,nv,-2)):
  if r["status"]!="pass_relaxed_mp2_torsion_scan_point" or v["status"] not in {"pass_independent_mp2_torsion_scan_point","pass_independent_relaxed_mp2_torsion_scan_point"} or r["rotor_id"]!=a.rotor_id or int(r["signed_step_index"])!=s:raise ValueError("terminal point mismatch")
  if v["scan_report"]["sha256"]!=sha(paths["positive_report" if s>0 else "negative_report"]):raise ValueError("validation binding mismatch")
 de={"positive_30":(float(pr["final_energy_hartree"])-a.frame_energy_hartree)*K,"negative_30":(float(nr["final_energy_hartree"])-a.frame_energy_hartree)*K}
 if not all(math.isfinite(x) and x<=10.0 for x in de.values()):raise ValueError("low-energy extension trigger absent")
 out=Path(a.output_dir).resolve();out.mkdir(parents=True,exist_ok=False)
 report={"schema_version":"1.0","report_type":"o6u_tier_a_low_energy_extension_authorization","status":"pass_tier_a_plus_minus45_extension_authorized","created_at_utc":datetime.now(timezone.utc).isoformat(),"rotor_id":a.rotor_id,"inputs":{k:rec(v) for k,v in paths.items()},"frame_energy_hartree":a.frame_energy_hartree,"terminal_relative_energies_kcal_mol":de,"trigger":"Both independently validated +/-30 terminals remain in the approximately 10 kcal/mol ffTK-relevant low-energy region.","authorized_signed_step_indices":[3,-3],"authorized_target_displacements_deg":[45,-45],"farther_extension_authorized":False,"chemistry_or_convergence_changed":False,"production_approved":False,"methodology_anchors":["10.1002/jcc.21367","10.1021/acs.jctc.5c00046","10.1002/jcc.23422","10.1063/5.0196848"],"disease_and_compound_anchors":["10.1136/ard.2010.134742","10.1002/art.30254","10.1158/1535-7163.MCT-10-0034","10.1038/s41594-024-01439-8"]}
 f=out/"O6U_TIER_A_LOW_ENERGY_EXTENSION_AUTHORIZATION.json";f.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n");print(json.dumps({"status":report["status"],"report":str(f),"sha256":sha(f),"relative_energies":de},sort_keys=True))
if __name__=="__main__":main()
