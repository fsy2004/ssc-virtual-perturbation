# O6U 配体参数化 全面审计与定案（2026-08-23 重建）

> 目的：解决项目内部"本地审计过时/互相矛盾"的问题。用户明确要求：
> "FFParam/QM 本地肯定有记录，现在项目太长了很多审计都过时并且矛盾，
> 需要把本地和服务器所有文件审计一遍，重新整理一下，然后更新记忆。"
> 本文档整合本地 + 服务器 `ssc_md_work` + release 三处的全部参数化证据，
> 逐条标注"已过时 / 已矛盾 / 仍生效"，并给出唯一终态。

---

## 0. 一句话定案（权威结论）

**O6U 参数化做了大量工作，但从未走到"最终批准"。**

- **采用进生产**的配体参数 = `common/toppar/O6U.itp`（CHARMM-GUI FF-Converter 2016-08-16 19:21 输出）——
  电荷经 FF-Converter **归一化**（非全零、总量=0、CGenFF 原子类型匹配），键/角/二面角力常数取自 **CGenFF/CHARMM36 官方 `forcefield.itp`**。
- **未完成**：FFParam/QM 电荷拟合与二面角扫描尝试了一整天（08-14~08-16），但**候选均未改善目标函数**（baseline objective 153.93，所有候选 >177），且**未生成 `ligand_parameter_record` 最终记录**、**无任何 `production_approved=true`**、16 日 12:00 后拟合目录无活动。
- **因此**：CGenFF 初始赋值的 penalty（max param 35.5 / max charge 32.824）是**已记录局限**，未被 FFParam/QM 拟合消除；但**不是"全零电荷"**（批注3那条错判不成立）。

---

## 1. 本地文件状态（注意：部分是早期/过时快照）

### 1.1 `inputs/ligand_parameterization/`
| 文件 | 状态字段 | 判定 |
|---|---|---|
| `O6U_PREPARAMETERIZATION_AUDIT.json` | `local_identity_pass_parameterization_blocked`；`md_parameterization_approved: false` | **仍生效**（身份过，参数化未批） |
| `O6U_CGENFF_INITIAL_ASSIGNMENT_AUDIT.json` | `pass_initial_assignment_only_requires_targeted_qm_validation`；`production_approved: false` | **仍生效**（初始赋值仅做身份；max penalty 35.5/32.824） |
| `O6U_CGENFF_MOL2_PREPARATION_AUDIT.json` | `pass_for_initial_cgenff_submission_only` | 仍生效（仅提交用） |
| `O6U_FFPARAM_TARGET_PLAN.json` | `draft_no_qm_target_generated`；`production_approved: false` | **过时**（08-10 状态；之后服务器上实际做了 QM 拟合目标生成，见 §2） |
| `O6U_SERVER_PARAMETERIZATION_HANDOFF.md` | 指令性文档 | 仍生效（描述"必须怎么做"，非"已完成"） |
| `O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.tsv` | 76 行映射 | 仍生效 |
| `O6U_neutral_*` | 输入 | 仍生效 |

**矛盾点**：`O6U_FFPARAM_TARGET_PLAN.json` 标 `draft_no_qm_target_generated`，但服务器 `charge_fit_*` / `bonded_torsion_*` 目录证明 QM 目标与拟合**实际已执行**（08-14~08-16）。本地这份是 08-10 早期的"尚未开始"状态，未随服务器更新。

### 1.2 项目根协议/门控/状态文档
| 文件 | 关键声明 | 判定 |
|---|---|---|
| `SERVER_EXECUTION_STATUS_20260810.md` | "Historical checkpoint"；"NO-GO ... until the O6U parameter model passes"; "No project production trajectory has been started" | **过时/矛盾**（08-10 说未跑生产，但 3×500ns 已跑完。唯一仍可参考的是工具/环境记录） |
| `LITERATURE_AND_FAILURE_GATES.md` | line20 "unvalidated CGenFF penalties 35.5/32.824 cannot enter the rebuild"; line62 "Ligand parameters ... Block all membrane/MD" | **部分过时**（事实"初始 penalty 未通过拟合消除"仍对；但"未 enter rebuild"与实际"已 enter MD"矛盾——实际上 MD 用了它） |
| `SERVER_READINESS.md` | "Current status: NO-GO"；line13-14 参数化门控项全未勾选 | **过时/矛盾**（NO-GO 检查清单从未更新为已跑；生产实际发生） |
| `LIGAND_PARAMETERIZATION_PROTOCOL.md` | fail-closed 验收门（是什么标准） | **仍生效**（作为"标准"而非"已达成"） |
| `FFPARAM_ACCESS_AND_VERSION_AUDIT_20260810.md` | "No ligand parameter is approved merely because FFParam installed"; FFParam 1.2.0 ≠ v2 | 仍生效 |
| `PROTONATION_MODEL.md` | 固定 O6U 中性/Asp257(-)/Asp385(0)；"predeclared fixed-protonation model" | 仍生效 |
| `MMGBSA_PBSA_SECONDARY_ANALYSIS_SENSITIVITY_WITHDRAWAL_20260822.md` | PB/GB 撤回 | 仍生效 |
| `MMGBSA_PBSA_SECONDARY_ANALYSIS_AMENDMENT_20260818.md` | 早期 PB/GB 修订 | 部分过时（已被 08-22 撤回覆盖） |

---

## 2. 服务器 `ssc_md_work/.../server_runs/o6u_parameterization/`（权威、最全）

**214 个子目录**，涵盖（按时间 08-10 → 08-16）：
- **MP2/水探针/QM 目标生成**：`ffparam_water_input_generation`、`water_*`、`mp2_*`（MP2 优化多次启动，卡在未收敛/需恢复）。
- **电荷拟合**：`charge_fit_*`（`target_universe` → `seeded_candidate_batch` → `seeded_refinement_batch` v2-v4 → `expanded_polar_site_*` → `whole_vector_*`）——**08-14 密集迭代一整天**，每个都有 `INDEPENDENT_VALIDATION`。
- **二面角/键参数自适应扫描**：`bonded_torsion_compound_specific_adaptive_scope`、`adaptive_scan_initial_cgenff_fitting_inputs` v1/v2 + consolidated。

### 关键终态证据
| 证据 | 内容 | 结论 |
|---|---|---|
| `adaptive_scan_fitting_manifest_20260816_v1/O6U_ADAPTIVE_SCAN_FITTING_MANIFEST.json` | 08-16 09:54 生成 | 二面角扫描 manifest（计划），非批准 |
| `charge_fit_expanded_polar_site_candidate_batch_frame0342_20260814_v1/O6U_CHARGE_FIT_SEEDED_BATCH_STATE.json` | `baseline_objective=153.93`；24 候选全部 `improved:false`，objective 177-188；`unfavourable_holdouts_non_attractive:true` | **电荷拟合候选未改善目标函数** |
| `find ... -iname "*ligand_parameter_record*" | 无（仅 `templates/ligand_parameter_record.template.json`） | **未生成最终参数记录** |
| `grep -rl '"production_approved":true'` | 全局 0 命中 | **无任何批准** |
| 拟合目录 `-newermt '2026-08-16 12:00'` | 无文件 | 拟合工作 08-16 中午后停止 |

---

## 3. 服务器 release（实际 MD 采用）

### 3.1 `common/toppar/O6U.itp`（24360 B，FF-Converter 输出，08-16 19:21）
- **无 `[parameters]` 段**——只有 `[atoms][bonds][pairs][angles][dihedrals][position_restraints]`，力常数引用 `forcefield.itp`。
- **电荷**：76 原子，总量 = 0.000（中性），非全零。与本地初始 rtf **不同**（见 §3.2），说明电荷经 FF-Converter 归一化/重分配。

### 3.2 `common/topol.top`：O6U 电荷对照（MD 实参 vs CGenFF 初始）
| 原子 | 初始 CGenFF (o6u__o6u.rtf) | **MD O6U.itp / forcefield.itp** |
|---|---|---|
| C8 (CG2O1) | +0.448 | **+0.333** |
| C23 (CG2R51) | +0.288 | **+0.326** |
| C24 (CG2R53) | +0.167 | **+0.298** |
| N2 (NG2R51) | −0.021 | **−0.140** |
| N4 (NG2S1) | −0.308 | **−0.330** |
| N5 (NG2R50) | −0.735 | **−0.601** |
| O (OG2D1) | −0.449 | **−0.481** |

`forcefield.itp`（CHARMM FF in GROMACS format）原子类型定义里 O6U 特有电荷 =
`CG2O1 0.333` / `CG2R51 -0.086` / `CG2R53 0.298` / `CG2R61 0.033` / `NG2R50 -0.601` / `NG2S1 -0.330` —— **与 O6U.itp 完全一致**，确认是一套。

### 3.3 键参数来源
`forcefield.itp` 头部 `CHARMM FF in GROMACS format` = **CGenFF/CHARMM36 官方域**。bonds 力常数（`CG2O1-NG2S1 3.096e5`、`CG2R51-NG2R50 3.347e5` 等）是 CGenFF 官方值，**未发现独立 QM 拟合替换的键参数**。

---

## 4. 矛盾/过时核心，及为何发生

1. **本地 ≠ 服务器**：本地 `server_runs/o6u_parameterization/` 只同步了 **water probe/MP2 前置子集**；服务器有完整 **charge_fit / bonded_torsion / adaptive_scan**（08-14~08-16）。本地审计文档停在 08-09~08-11，未反映后续拟合工作。
2. **"未做" vs "做了但不通过"**：本地 `O6U_FFPARAM_TARGET_PLAN.json` 标 `draft_no_qm_target_generated`（=从未开始 QM）；但服务器证明 QM 拟合**实际执行了但候选未改善目标**。两条说法都"对"但含义不同——**本地应更新为"QM 拟合已尝试，未收敛到改善，未批准"**。
3. **"NO-GO 未跑生产" vs "已跑 3×500ns"**：`SERVER_EXECUTION_STATUS_20260810` / `SERVER_READINESS` 是**生产前**的 NO-GO 检查点；生产实际发生（08-16 后），但**这些检查清单从未被更新为"已跑"**，造成"文档说 NO-GO、实际跑了"的观感矛盾。

---

## 5. 终态声明（用于稿件+记忆）

> O6U（nirogacestat）使用 CHARMM-GUI/CGenFF 复合物特定参数：电荷经 FF-Converter 归一化（总量 0，非全零，CGenFF 原子类型匹配），键/角/二面角力常数源自 CGenFF/CHARMM36 官方域（`forcefield.itp`）。初始 CGenFF 赋值 penalty（max param 35.5 / max charge 32.824，超出常用 ~10 建议验证阈值）为**已记录局限**——FFParam/QM 电荷拟合与二面角扫描已做但**候选未能改善目标函数（baseline 153.9 vs 候选 ≥177），未生成最终参数记录、未获生产批准（无 production_approved=true）**。因此 MD 结果定位为**在采用参数下的探索性结构稳定性观察**，不支撑结合/亲和/药效/occupancy 结论；不声称 FFParam/QM 全面优化并最终批准。质子化状态为 predeclared fixed-protonation model（O6U 中性 / PSEN1 Asp257- / Asp385 0）。

---

## 6. 需要后续处理（供用户确认）
- [ ] 更新 `O6U_FFPARAM_TARGET_PLAN.json` 状态为"QM 拟合已尝试、未改善目标、未批准"（消除 draft_no_qm_target_generated 的过时）。
- [ ] 是否在 `LITERATURE_AND_FAILURE_GATES.md` / `SERVER_READINESS.md` 补一行"生产已于 08-16 执行，参数化未最终批准（见 §5）"以消除"NO-GO 但跑了"的矛盾。
- [ ] 稿件 Methods C7 配体描述按 §5 口径写（已记录局限，非 FFParam/QM 完全优化）。
