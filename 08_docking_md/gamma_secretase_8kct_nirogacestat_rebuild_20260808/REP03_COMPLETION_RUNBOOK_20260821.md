# rep03 完成后执行 runbook（2026-08-21 冻结，rep03 完成时按序执行）

> 每步在服务器 release 目录（`/root/autodl-tmp/o6u_md_release_3x500ns_v4`）执行。
> 分析 Python 环境：`/root/autodl-tmp/envs/ssc_md_analysis_py311/bin/python`
> GROMACS：`/root/GROMACS-2025.2/bin/gmx`
> 所有输出目录必须新建；raw 文件永不修改。

## Step 0：确认状态

```bash
python scripts/status_new_md_run.py   # 本地运行；确认 rep03 FINISHED=yes、无 gmx 进程、哈希匹配
```

## Step 1：rep03 终点门控

```bash
python scripts/validate_start_next_production.py --completed rep03
# 生成 rep03/PRODUCTION_COMPLETION_500NS.json（status=pass 必须）
```

## Step 2：三重复 raw 完整性验证

```bash
python scripts/validate_md_outputs.py --manifest config/study_manifest.json --phase production --strict --report reports/production_output_validation.json
# 需要 study_manifest.json 中 realizations 哈希已填（rep03 完成后填充）
```

## Step 3：手动 PBC 预演（生成 05_centered 供分析）

在 release 目录下按 make_analysis_trajectories 计划的命令手工执行（whole→cluster→nojump→center→fit），输出到临时目录 `analysis/prep_<rep>/`；每重复处理完删除中间 01/02/04/05，保留 06_fitted/07_fixed/命令记录/日志/哈希。

## Step 4：分析脚本（结构/膜/能量）

```bash
PY=ssc_md_analysis_py311/bin/python
$PY scripts/analyze_primary_structure_mdanalysis.py --manifest config/primary_postprocessing_manifest.json --output-root analysis_primary/
$PY scripts/analyze_membrane_qc_mdanalysis.py --manifest config/primary_postprocessing_manifest.json --output-root analysis_primary/
$PY scripts/gmx_energy_qc.py --manifest config/primary_postprocessing_manifest.json --output-root analysis_primary/ --mode extract
# 前置：primary_postprocessing_manifest.json 的 realizations 哈希填完 + approval_status=approved
# 预期：membrane COMPLETE preproduction_status=blocked_external_membrane_metrics（APL/gorder NO-GO）
```

## Step 5：QC/stationarity 评估与密封

1. 人工审查结构/膜/能量 COMPLETE.json 与 block 数据（对照 ANALYSIS_CONFIG_FREEZE_20260821.md 阈值）
2. 填写 qc_draft（criterion evidence：帧数/5 块/4 门控）
3. `python scripts/validate_qc_stationarity_report.py --manifest config/study_manifest.json --analysis-plan config/analysis_plan.json --raw-output-report reports/production_output_validation.json --report <qc_draft> --seal-output reports/qc_sealed.json`
4. 把 qc_sealed.json 的 SHA-256 填入 analysis_plan.json 的 eligibility_gate.qc_and_stationarity_report_sha256

## Step 6：正式轨迹链

```bash
python scripts/make_analysis_trajectories.py --manifest config/study_manifest.json --analysis-plan config/analysis_plan.json --production-qc-report reports/qc_sealed.json --gmx /root/GROMACS-2025.2/bin/gmx --execute
# 生成 analysis/trajectories/... 01-07 + mindist + PBC 距离不变性报告；每重复完成后清理 01/02/04/05 中间文件（磁盘约束）
```

## Step 7：独立验证

```bash
$PY scripts/validate_primary_postprocessing.py --manifest config/primary_postprocessing_manifest.json --output-root analysis_primary/ [--review-dispositions <spike处置>]
# 预期 claim_gate：若外部膜指标未验证 → blocked_inconclusive（如实报告证据缺口）
```

## Step 8：PB/GB 次要分析（本节点 CPU 模式）

1. `seal_secondary_endpoint_all_three_gate.py`（全部门控密封）
2. `install_endpoint_preprocess_env.sh <prefix>`（如已存在 ssc_md_analysis_py311 则复用/核对）
3. 每重复 `prepare_secondary_endpoint_energy_inputs.py`（300 帧 + 3 帧 canary）
4. `freeze_endpoint_energy_membrane_geometry.py`（mthick 冻结；需膜 QC preproduction_status=pass——**若 blocked 则按协议 No-Go 报告**）
5. `build_endpoint_energy_cpu_migration_package.py`（迁移包）
6. `install_gmx_mmpbsa_1_6_5_cpu.sh <prefix>`（**安装 gmx_MMPBSA 1.6.5 属既定方案，执行前向用户报告**）
7. `capture_gmx_mmpbsa_toolchain.py`、`run_gmx_mmpbsa_canary.py --execute`（3 帧 canary）
8. `collect_endpoint_cpu_inventory.py`、`plan_secondary_endpoint_resources.py`
9. `run_secondary_endpoint_energy_cpu.py`（4 模型 × 3 重复）
10. `summarize_secondary_endpoint_energy.py`（描述性汇总，无 p 值/排序/亲和力声明）

## Step 9：后续（主图/稿件/交付）

见交付计划；全部数字回查后写稿件。

## 关键约束提醒

- 阈值已预注册（ANALYSIS_CONFIG_FREEZE_20260821.md），不得结果导向修改
- 外部膜指标未验证时：claim gate=blocked_inconclusive；PB/GB 的 mthick 冻结要求膜 QC pass——若 blocked，PB/GB 的几何冻结步骤按协议 No-Go（或如实报告后评估）
- CUDA 恢复证据（rep01/02/03 audit）必须纳入连续性分析
- 不覆盖 raw；不删除帧；不改种子/协议
