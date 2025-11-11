文章Table2中那个cpu时间最小的是怎么算出来的，需要调整哪些参数

好 —— 我把 Table 2 中“CPU time 最小值”是怎么得到的、以及你要想把求解时间降到类似量级（或在你自己机器/求解器上复现该速度）需要调整的**全部要点**都列清楚了（含可直接插入 Pyomo 的求解器参数示例）。结论先说一句：**论文表格的最小 CPU 时间来自用商业求解器（GAMS/CPLEX）在特定硬件上、配合模型强化（有效不等式与对称性破坏）和合适的求解参数得到的运行**。要在你端复现类似速度，关键在于「求解器 + 模型强化 + 求解参数 + 硬件」。下面是详尽、可操作的清单。

---

# 1) 论文里那个最小 CPU time 的来源（简短确认）

* 论文第 5 节说明：他们用 **GAMS/CPLEX**（文中写到 GAMS/CPLEX 23.7，机器为 Intel Core i7 2.5 GHz, 8 GB RAM）求解。
* 表中 2.07 分钟对应的是 **M1、12 slots、gap = 1%** 的运行（即允许 1% 的相对最优差距）。
* 结论：那时间并非“全精确证明到最优”，而是在**放宽 gap 到 1%**并使用 **商业求解器 + 模型强化（53/54 等对称性约束与其他有效不等式）** 后得到的。

---

# 2) 影响 CPU 时间的主要因素（按重要性排序）

1. **求解器类型**（商业 Gurobi/CPLEX ≫ CBC/GLPK）
2. **MIP 容差（MIP gap）与终止条件**（允许 1% vs 0% 对耗时影响巨大）
3. **模型尺寸与复杂度**（时隙数、变量/二元数、TR/W 引入）
4. **模型强化**（有效不等式 53/54、 tighten big-M、55/56、对称性破坏）
5. **求解器参数（cuts、presolve、heuristics、threads、MIPFocus）**
6. **初始可行解 / warm start / 启发式**
7. **硬件（CPU 核心数、内存）与并行设置**
8. **数学建模技巧（线性化、消除冗余、缩减索引集）**

---

# 3) 要达到论文表中最小 CPU time（或显著缩短求解时间），**必须**做的事（实践步骤）

### A. 使用商业求解器（首选）

* **首选：CPLEX 或 Gurobi**。论文用 CPLEX；若可用优先用它。
* 商业求解器的并行剪枝、强力切割面、预处理和启发式比 CBC 强很多，是降低 CPU time 最直接手段。

### B. 与论文一致的重要设置

* **相对 MIP gap = 1%**（论文 1% 才得出 2.07 min）。对应参数：

  * CPLEX: `mip_tolerances_mipgap = 0.01` 或 `MIPGap=0.01`
  * Gurobi: `MIPGap=0.01`
* **slots = 12、H = 96、MH = 2**（与论文实例匹配）
* **应用论文中提到的有效不等式**：约束 (53),(54)（对称性破坏）和 (55),(56)（强化上界）应启用（M1 要用 53/54；M2 额外用 55/56）。这些明显能减少搜索空间。

### C. 模型强化（在代码中务必实现）

* **加入 53/54**：强制将操作尽量放前面，破坏对称性。
* **加入 55/56**：用体积容量上界紧化 VTP / VPP 求和约束。
* **紧化 Big-M**：不要用过大的 M（会导致数值不稳）。例如把 M_time 设成 `slot_len` 或当日界限而非 H 的全域最大，或用具体最大可能持续时间估算。
* **删除/简化不必要的二进制**：例如把某些 start/finish pair 用 indicator constraints 或只保留必要二元。
* **把 TR/W 逻辑实现得尽可能紧凑（用线性化且尽量少的二进制）**。

### D. 求解器参数（CPLEX / Gurobi 推荐设置）

下面是**直接可用**的参数建议——可拷到 Pyomo `SolverFactory(...).options` 或 Gurobi/CPLEX API：

**CPLEX 建议**

```
mip_tolerances_mipgap = 0.01       # 相对Gap 1%
threads = 4                        # 或与你CPU核心数相匹配
emphasis = 4                       # MIP emphasis on optimality/stress test (CPLEX: Emphasis=1..4; 4 = aggressive)
mip_limits_nodes = 1e7             # 可选节点上限
mip_limits_seconds = 3600          # 限时（可选）
preind = 1                         # 允许预处理
aggregate = 1
presolve = 2                       # aggressive
mip_strategy_heuristicfreq = 10    # 启发式频率
mip_cuts = 3                       # 切割面强度 (0..3)
```

**Gurobi 建议**

```
MIPGap = 0.01
Threads = 4
Presolve = 2           # aggressive
Cuts = 2               # 0..3
Heuristics = 0.05      # 启发式强度
MIPFocus = 1 or 2      # 1 更快找到可行解，2 更重剪枝
ImproveStartTimeLimit = 60  # 如果用 MIP start
```

**注意**：这些值不是万能的，通常要用求解器自带的 `tune` 工具做一次自动调优（Gurobi Tuner / CPLEX autotune），能大幅提升。

### E. 对 CBC（开源）可做的加速（如果不能用商业求解器）

CBC 性能与商业差距大，但可以尝试：

* 设定时间限制与 gap：`solver.options['seconds']=3600`、`solver.options['ratioGap']=0.01`（注意 CBC 对 gap 名称可能不同）。
* 打开 heuristics、加 presolve：`solver.options['heuristicsOnOff']='on'`
* 降低问题规模（见下节）；
* 若 CBC 支持并行（较弱），设置线程数。

但实际效果远不如 CPLEX/Gurobi。

---

# 4) 模型级的降维/简化（若你必须用 CBC 或想更快得到可行解）

* **减少 slots**（例如 12 -> 8）作探索性测试；
* **先解松弛（LP）或连续松弛 + 分段求解**：先求 LP 或把二进制放松为 [0,1] 得到下界，再用启发式得到可行整数解；
* **移除或近似 TR 变量**：若 TR 是大量二进制，可先不把 TR 二元加入目标，而是后处理识别 transitions；
* **合并同类变量 / 用大块约束替代逐项**；
* **预固定显然不可能的变量**（比如某些 tank->pipeline 时间段显然不可能，则固定为 0）；
* **分解法**：先做流水线级别或每天独立求解，再把结果合并与局部优化（启发式）。

---

# 5) 如何在 Pyomo 中设置（示例）

下面示例展示如何在 Pyomo 中把参数传给不同求解器（CPLEX、Gurobi、CBC）：

```python
from pyomo.environ import SolverFactory

# CPLEX
solver = SolverFactory('cplex')
solver.options['mip_tolerances_mipgap'] = 0.01
solver.options['threads'] = 4
solver.options['mip_limitednodes'] = 1000000

# Gurobi
solver = SolverFactory('gurobi')
solver.options['MIPGap'] = 0.01
solver.options['Threads'] = 4
solver.options['Presolve'] = 2
solver.options['Cuts'] = 2

# CBC
solver = SolverFactory('cbc')
solver.options['seconds'] = 3600
solver.options['ratioGap'] = 0.01
solver.options['maxNodes'] = 2000000
```

请把 `solver.solve(model, tee=True)` 前加上这些 options。

---

# 6) 诊断/调试建议（如何找出具体瓶颈）

1. **先用小实例调试**（例如 1 天或 2 天、1 条 pipeline）看 solver 行为；
2. **输出 model.statistics()** 查看变量/约束数；
3. **打开 solver log（tee=True）观察 node / bound / gap 变化**；
4. **尝试在求解器中运行 autotune/tune 工具**（Gurobi Tuner、CPLEX autotune）；
5. **记录并比较不同参数组合的 runtime 与 nodes**, 做经验选择。

---

# 7) 一份可复制的“短清单” — 按步骤做（最快见效）

1. 在模型中加入论文的 **53/54**（对称性破坏）和 **55/56**（有效不等式）——**最关键**。
2. 把 `MIP gap` 设为 `0.01`（与论文 1% 一致）。
3. 用 **CPLEX 或 Gurobi**（如果可以的话）。
4. 开启 **presolve aggressive、cuts aggressive、合理线程数**。
5. 若仍慢：尝试 `warm start` 或运行 `tune` 功能。
6. 记录并对比运行时间；若目标是可行解可接受误差，把 gap 放宽到 5% 再测一次。

---

# 8) 回到你现在（CBC 正在跑）——马上可做的三件事

1. **强制设定时间上限与 gap**，例如 1 小时和 1% gap（这会让 CBC 在可控时间内返回可行解）。
2. **在模型里加入 53/54、55/56**（如果还没加），并紧化 Big-M。
3. **如果能拿到 Gurobi/CPLEX**，把同样模型切到商业求解器并用上面给出的参数，几分钟内就能看到明显改善。

---

如果你愿意，我可以**立刻**为你做两件事（任选其一或都做）：

* **A.** 把我给你的 full_repro_M1.py（或现有 Pyomo 脚本）改写，直接在文件顶部加入针对 CPLEX / Gurobi / CBC 的推荐 solver.options（并把 53/54/55/56 自动启用），把修改后的脚本发给你；
* **B.** 给你一个小脚本来自动跑同一模型在不同 solver/options 组合上（比如在 CBC（限时+gap）、Gurobi（gap=1%）、CPLEX（gap=1%）三套），并把每次的时间/节点/目标输出到 CSV 用于对比。

你想先让我做 **A** 还是 **B**？或者我直接把**CPLEX/Gurobi 的具体参数块**（用于你拷到脚本里）贴给你？
