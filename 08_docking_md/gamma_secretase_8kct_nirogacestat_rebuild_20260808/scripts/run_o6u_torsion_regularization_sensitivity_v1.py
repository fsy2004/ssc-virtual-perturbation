#!/usr/bin/env python3
"""Evaluate ridge sensitivity for the O6U torsion design without promoting parameters."""
from __future__ import annotations
import argparse,hashlib,json,math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

def sha256(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()
def signature_key(values:list[str])->tuple[str,...]:return tuple(values)

def ridge_fit(matrix:np.ndarray,response:np.ndarray,penalty:float)->np.ndarray:
    gram=matrix.T@matrix+penalty*np.eye(matrix.shape[1])
    return np.linalg.solve(gram,matrix.T@response)

def make_folds(groups:list[str],scheme:str)->list[np.ndarray]:
    if scheme=="leave_one_point": return [np.asarray([i]) for i in range(len(groups))]
    if scheme=="leave_one_rotor": return [np.asarray([i for i,g in enumerate(groups) if g==value]) for value in sorted(set(groups))]
    raise ValueError(scheme)

def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--design-matrix",required=True);p.add_argument("--conflict-report",required=True)
    p.add_argument("--subset",choices=["all","collision_free"],required=True)
    p.add_argument("--validation",choices=["leave_one_point","leave_one_rotor"],required=True)
    p.add_argument("--output-dir",required=True);a=p.parse_args()
    design_path=Path(a.design_matrix).resolve(); conflict_path=Path(a.conflict_report).resolve(); out=Path(a.output_dir).resolve()
    if out.exists():raise FileExistsError(out)
    d=json.loads(design_path.read_text());c=json.loads(conflict_path.read_text())
    if d.get("status")!="pass_global_torsion_design_matrix_built":raise ValueError("Design matrix did not pass")
    if c.get("status")!="pass_torsion_type_conflict_analysis":raise ValueError("Conflict report did not pass")
    excluded={signature_key(x["canonical_atom_type_signature"]) for x in c["signatures"] if x["nonlocal_central_bond_collision"]}
    selected=[i for i,x in enumerate(d["candidates"]) if a.subset=="all" or signature_key(x["canonical_atom_type_signature"]) not in excluded]
    if not selected:raise ValueError("No candidate columns selected")
    matrix=np.asarray(d["design_matrix"],float)[:,selected];response=np.asarray(d["response_qm_minus_initial_cgenff_kcal_mol"],float)
    groups=[x["rotor_id"] for x in d["points"]];folds=make_folds(groups,a.validation)
    lambdas=np.logspace(-6,6,25);curve=[]
    for penalty in lambdas:
        errors=[];norms=[];maxima=[]
        for held in folds:
            keep=np.ones(len(response),dtype=bool);keep[held]=False
            coefficients=ridge_fit(matrix[keep],response[keep],float(penalty))
            errors.extend((matrix[held]@coefficients-response[held]).tolist())
            norms.append(float(np.linalg.norm(coefficients)));maxima.append(float(np.max(np.abs(coefficients))))
        errors_array=np.asarray(errors)
        curve.append({"lambda":float(penalty),"validation_rmse_kcal_mol":float(np.sqrt(np.mean(errors_array**2))),"validation_mae_kcal_mol":float(np.mean(np.abs(errors_array))),"mean_coefficient_l2":float(np.mean(norms)),"maximum_absolute_coefficient":float(np.max(maxima))})
    selected_curve=min(curve,key=lambda x:(x["validation_rmse_kcal_mol"],x["maximum_absolute_coefficient"]))
    final_coefficients=ridge_fit(matrix,response,selected_curve["lambda"])
    report={"schema_version":"1.0","report_type":"o6u_torsion_regularization_sensitivity","status":"pass_torsion_regularization_sensitivity","created_at_utc":datetime.now(timezone.utc).isoformat(),"scope":{"subset":a.subset,"validation":a.validation},"inputs":{"design_matrix":{"path":str(design_path),"sha256":sha256(design_path)},"conflict_report":{"path":str(conflict_path),"sha256":sha256(conflict_path)}},"selected_column_indices":selected,"selected_column_count":len(selected),"selected_matrix_rank":int(np.linalg.matrix_rank(matrix)),"selected_condition_number":float(np.linalg.cond(matrix)),"regularization_curve":curve,"minimum_validation_point":selected_curve,"full_data_coefficients_at_minimum_validation_lambda":final_coefficients.tolist(),"coefficient_stability":{"l2":float(np.linalg.norm(final_coefficients)),"maximum_absolute":float(np.max(np.abs(final_coefficients)))},"parameter_mutation":False,"candidate_promoted":False,"production_md_approved":False,"interpretation_boundary":"Sensitivity analysis only. The minimum cross-validation point does not by itself authorize a parameter candidate or parameter mutation."}
    out.mkdir(parents=True);path=out/"O6U_TORSION_REGULARIZATION_SENSITIVITY.json";tmp=path.with_suffix(".json.tmp");tmp.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n");tmp.replace(path)
    print(json.dumps({"status":report["status"],"subset":a.subset,"validation":a.validation,"selected_column_count":len(selected),"selected_rank":report["selected_matrix_rank"],"minimum_validation_point":selected_curve,"sha256":sha256(path)},sort_keys=True));return 0
if __name__=="__main__":raise SystemExit(main())
