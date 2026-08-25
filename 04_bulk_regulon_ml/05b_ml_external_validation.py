#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""05b · ML 外验(现代化增量): regulon 活性特征 -> SSc vs HC 分类, 留一队列外验 AUROC。
★铁律([[reference_dl_ai_strategy]]): 复杂/FM 法必配简单基线, 打不过基线=白搭 -> 同时跑
TabPFN(2024 表格 FM, 打得过 GBDT) + Logistic 基线, 报两者 AUROC, 让审稿人看谁赢。
可解释: permutation importance(哪些 regulon/lead 驱动分类) -> 佐证 HES1/FOSB。
支撑层非头条(association, 非因果)。读 05 缓存的 {GSE}_regulon_activity.csv + _meta.csv。
env: npc (+ optional tabpfn)。用法: SSC_OUT_05=out_05_bulk python 05b_ml_external_validation.py"""
import os, sys, glob, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from sklearn.inspection import permutation_importance

OUT05 = os.environ.get("SSC_OUT_05", "out_05_bulk")
FEAT  = os.path.join(OUT05, "features")
OUT   = os.path.join(OUT05, "ml"); os.makedirs(OUT, exist_ok=True)
LEADS = ["HES1", "FOSB", "SMAD3", "MEF2C"]
SEED  = 0
def log(m): print(f"[05b-ml] {m}", flush=True)

# ---- 1. 汇集各队列 regulon 活性(样本×regulon) + SSc/HC 标签 ----
acts, metas = {}, {}
# ---- documented cohort QC exclusion (2026-07-10; not cherry-picking) ----
# GSE76885 = longitudinal MMF-treatment / registry cohort (172 SSc : 18 HC). Its case-control
# contrast is confounded: the 18 "HC" are registry normals with HIGH fibrotic signal (fibro_load
# 0.50 vs treated-SSc -0.05), so even SMAD3 (positive control) is mis-oriented -> it FAILS the
# positive-control-orientation QC (the same "pc_pass" gate the meta-analysis uses to drop it).
# It is a treatment-trajectory cohort, NOT a cross-sectional diagnostic cohort, so it does not
# belong in the SSc-vs-HC transfer benchmark. (Within its SSc patients the regulons still track
# fibrosis: SMAD3 rho=0.73 vs fibro_load -> the panel is valid; only the diagnostic contrast is
# confounded.) Analysed separately, not in LOCO.
EXCLUDE_DIAGNOSTIC = {"GSE76885"}
for f in sorted(glob.glob(os.path.join(FEAT, "*_regulon_activity.csv"))):
    gse = os.path.basename(f).split("_regulon")[0]
    if gse in EXCLUDE_DIAGNOSTIC:
        log(f"排除 {gse}(治疗/登记队列, 阳性对照方向 QC 不合格 -> 单独分析, 不进诊断 LOCO)"); continue
    a = pd.read_csv(f, index_col=0)
    mf = os.path.join(FEAT, f"{gse}_meta.csv")
    if not os.path.exists(mf):
        log(f"跳过 {gse}(无 meta)"); continue
    m = pd.read_csv(mf, index_col=0)
    # 标签列: SSc/status/group 里找 0/1 或 SSc/HC
    lab_col = next((c for c in ["grp", "SSc", "status", "group", "condition", "disease"] if c in m.columns), None)
    if lab_col is None: log(f"跳过 {gse}(无标签列)"); continue
    y = m[lab_col].astype(str).str.upper().map(
        lambda s: 1 if ("SSC" in s or s in ("1", "TRUE", "CASE", "DISEASE")) else (0 if ("HC" in s or "CTRL" in s or "CONTROL" in s or s in ("0", "FALSE", "NORMAL", "HEALTHY")) else np.nan))
    common = a.index.intersection(y.dropna().index)
    if len(common) < 8: log(f"跳过 {gse}(n={len(common)}<8)"); continue
    acts[gse] = a.loc[common]; metas[gse] = y.loc[common].astype(int)
    log(f"{gse}: {len(common)} 样本 ({int(metas[gse].sum())} SSc / {int((metas[gse]==0).sum())} HC), {a.shape[1]} regulon")
if len(acts) < 2:
    log("外验需 >=2 队列; 现有 {} -> 退出(等 05 产出更多队列特征)".format(len(acts))); sys.exit(0)

# 共有 regulon 特征
feats = sorted(set.intersection(*[set(a.columns) for a in acts.values()]))
log(f"共有 regulon 特征: {len(feats)}")
cohorts = [c for c in acts if metas[c].nunique() == 2]
log(f"诊断 LOCO 队列(2 类平衡): {cohorts}")

# ---- 2. 留一队列外验(LOCO): 训 N-1 队列, 测留出队列 ----
def fit_predict(model_name, Xtr, ytr, Xte):
    if model_name == "logistic":
        sc = StandardScaler().fit(Xtr)
        clf = LogisticRegression(max_iter=2000, C=1.0, random_state=SEED).fit(sc.transform(Xtr), ytr)
        return clf.predict_proba(sc.transform(Xte))[:, 1], clf
    if model_name == "tabpfn":
        from tabpfn import TabPFNClassifier
        clf = TabPFNClassifier(device="cpu", ignore_pretraining_limits=(Xtr.shape[1] > 500))   # ★>500 特征: 抑制 TabPFN 预训练上限校验(ignore_pretraining_limits, 否则 fit 抛 TabPFNValidationError); 此时超出预训练域(OOD), 其 AUROC 作 OOD 参考、不当干净 FM 基线; 2026-07-08 修
        clf.fit(Xtr, ytr)
        return clf.predict_proba(Xte)[:, 1], clf
    raise ValueError(model_name)

rows = []
# ★2026-07-07 修: 原守卫只测 import(成功即入 models), 但 tabpfn 8.x 运行期需 license/权重,
# fit() 抛 TabPFNLicenseError -> 逐折静默失败, 汇总只剩基线却像跑过 TabPFN 对比。改为真跑
# 一次微型 fit 验证可用性, 不可用则干净排除(诚实, 不污染结论)。
has_tabpfn = False
try:
    from tabpfn import TabPFNClassifier as _TPC
    _t = _TPC(device="cpu")
    _t.fit(np.array([[0., 1.], [1., 0.], [0., 0.], [1., 1.]]), np.array([0, 1, 0, 1]))
    _t.predict_proba(np.array([[0.5, 0.5]]))
    has_tabpfn = True
except Exception as _e:
    log(f"★tabpfn 不可用({type(_e).__name__}: {_e}) -> 仅跑 Logistic 基线(如需 TabPFN 对比: 设 TABPFN_TOKEN 或预缓存权重)")
models = ["logistic"] + (["tabpfn"] if has_tabpfn else [])
for held in cohorts:
    tr = [c for c in cohorts if c != held]
    Xtr = pd.concat([acts[c][feats] for c in tr]).values
    ytr = pd.concat([metas[c] for c in tr]).values
    Xte = acts[held][feats].values; yte = metas[held].values
    if len(np.unique(yte)) < 2:  # 留出队列单类无法算 AUROC
        log(f"留出 {held} 单类, 跳过 AUROC"); continue
    for mdl in models:
        try:
            p, _ = fit_predict(mdl, Xtr, ytr, Xte)
            auc = roc_auc_score(yte, p)
            rows.append(dict(held_out=held, model=mdl, n_test=len(yte), auroc=round(auc, 3)))
            log(f"  留出 {held} · {mdl}: AUROC={auc:.3f}")
        except Exception as e:
            log(f"  留出 {held} · {mdl} 失败: {e}")
res = pd.DataFrame(rows); res.to_csv(os.path.join(OUT, "loco_auroc.csv"), index=False)
# 汇总: 各模型平均 AUROC(TabPFN 必须 >= 基线才有意义)
summ = res.groupby("model")["auroc"].agg(["mean", "std", "min"]).round(3)
summ.to_csv(os.path.join(OUT, "loco_auroc_summary.csv"))
log("留一队列外验汇总(★TabPFN 须 >= Logistic 基线):\n" + summ.to_string())

# ---- 3. 可解释: permutation importance(哪些 regulon 驱动分类) ----
# ★2026-07-07 修: 原在全数据(=拟合所用同一批)上算 PI = 重代入/记忆化(sklearn 文档: PI 须在留出集上算)。
# 改为逐 LOCO fold 在留出队列上算 PI, 跨 fold 聚合 -> 反映泛化而非记忆。
fold_imps = []
for held in cohorts:
    tr = [c for c in cohorts if c != held]
    Xtr = pd.concat([acts[c][feats] for c in tr]).values; ytr = pd.concat([metas[c] for c in tr]).values
    Xte = acts[held][feats].values; yte = metas[held].values
    if len(np.unique(yte)) < 2:            # 留出队列单类 -> 无法算 roc_auc PI, 跳过该 fold
        continue
    scl = StandardScaler().fit(Xtr)
    clf = LogisticRegression(max_iter=2000, random_state=SEED).fit(scl.transform(Xtr), ytr)
    p = permutation_importance(clf, scl.transform(Xte), yte, n_repeats=20, random_state=SEED, scoring="roc_auc")
    fold_imps.append(p.importances_mean)
if fold_imps:
    imp_mat = np.vstack(fold_imps)
    imp = pd.DataFrame({"regulon": feats, "importance": imp_mat.mean(0), "std": imp_mat.std(0)}) \
            .sort_values("importance", ascending=False)
else:                                       # 无可用留出 fold -> 空表, 不伪造重要性
    imp = pd.DataFrame({"regulon": feats, "importance": np.nan, "std": np.nan})
imp.to_csv(os.path.join(OUT, "regulon_importance.csv"), index=False)
log("top10 判别 regulon:\n" + imp.head(10).to_string(index=False))
log("★lead(HES1/FOSB) 排名: " + ", ".join(
    f"{g}=#{int(imp.reset_index(drop=True).index[imp['regulon']==g][0])+1}" for g in LEADS if g in imp['regulon'].values))
log(f"done -> {OUT}/ (loco_auroc, summary, regulon_importance)")
