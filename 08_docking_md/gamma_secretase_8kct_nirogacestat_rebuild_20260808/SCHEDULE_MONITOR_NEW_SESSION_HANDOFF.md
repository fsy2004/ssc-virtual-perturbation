# rep03 定时监控 · 新会话交接指令（dsh-schedule 版）

> 用途：您新建一个 DSH 会话（root agent）后，用它原生 `schedule_create` 每 30 分钟
> 检查并报告 rep03 进度。新会话创建时机须在 @deepseek-ai/dsh-schedule 插件加载之后
> （当前 profile 已挂载 schedule，故任何新开会话都会获得这三个工具）。

## 一、给新会话的完整任务指令（可直接粘贴)

```
你负责监控服务器上正在运行的分子动力学生产任务 rep03，直到其到达 500 ns 完成。
请你：

1. 调用 schedule_create 创建一个**固定间隔**提醒：every_seconds=1800，
   prompt="检查 rep03 MD 生产进度并报告：运行 python scripts/status_new_md_run.py，
   读取 rep03 的 ps 值、gmx PID、GPU、磁盘、OOM；每 30 分钟向用户报告一次状态；
   若 rep03 达到 500 ns（FINISHED=yes）或出现异常（OOM>0、gmx 进程消失、哈希不匹配），
   立即报告并停止该 schedule。"

2. 在每 30 分钟的提醒里执行上述检查，向用户报告一行摘要（时间/ps/gmx/gpu/disk/oom），
   并读取本地状态文件获取最新值（见下）。

3. 参考命令与路径：
   - 检查脚本: C:\Users\fsy\Desktop\SSc_project\06_code_reproducibility\08_docking_md\gamma_secretase_8kct_nirogacestat_rebuild_20260808\scripts\status_new_md_run.py
   - 实时状态文件: ...\reports\rep03_monitor\latest_status.txt
   - 完整日志: ...\reports\rep03_monitor\rep03_monitor.log
   - 服务器端点、哈希、监控器逻辑见该目录 scripts\monitor_rep03_to_completion.py

4. 当 rep03 完成（500 ns）时，按顺序执行终点门控，脚本为：
   python scripts\validate_start_next_production.py --completed rep03
   （后续分析链见 ...\REP03_COMPLETION_RUNBOOK_20260821.md）
```

## 二、当前状态（2026-08-22，交接时点）

- rep03：**460.32 ns / 500 ns（92.1%）**，gmx PID 210062 运行中
- GPU：RTX 5090，66%/76°C；磁盘 238G 剩余；cgroup OOM/oom_kill = 0
- 冻结哈希全部匹配（archive/manifest/canary/release）
- 端点：connect.westc.seetacloud.com:29621（凭据由 scripts/new_md_server.py 从用户密码文件读取，勿复制）
- 剩余约 40 ns ≈ 3.5 小时

## 三、说明

- 若新会话在 dsh-schedule 加载前创建，则看不到 schedule_create——请确保是新开
  （不是从当前会话 fork/subagent，子代理不是 root agent，不会有该工具）。
- 提醒是 session-local：新会话须保持 live，否则提醒延迟到恢复后补发。
- 若新会话无法创建/用 schedule，可退回常驻监控+心跳链方案（当前会话已在运行）。
