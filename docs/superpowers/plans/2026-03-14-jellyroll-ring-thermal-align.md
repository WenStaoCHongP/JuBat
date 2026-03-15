# Jellyroll vs 圆环热模型尺度对齐修改计划

> **For agentic workers:** REQUIRED: Use superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在验证脚本中让 Ring2D 热模型与 Jellyroll 使用同一套无量纲尺度（scale.L、scale.q、scale.h 等），并修正热源桥接方式，使两条温度曲线形态一致（均单调上升、量级接近），消除圆环“上冲后下降”的差异。

**Architecture:** 在 `run_ring_simulation` 中，在调用 `SetCase(param_ring, opt)` 之前，将 Jellyroll 的 `param_dim.scale` 中与热求解相关的字段复制到 `param_ring.scale`，使 Ring 与 Jellyroll 使用相同参考长度、热源尺度与 Biot 数；热源输入改为直接使用体热源 [W/m³]，不再乘以 `cell.area`。

**Tech Stack:** Julia, JuBat (src/SetParams.jl, example/热模块验证/jellyroll_vs_ring_thermal_compare.jl)

**参考:** `md/05_热模型_二维分布式.md`（边界条件：`cool_method="none"` 仅表示不添加 z 向冷却，侧面仍施加对流）；`src/SetParams.jl` 中 `scale.L = PE.thickness + NE.thickness + SP.thickness`，`scale.q = P_ref/L^3`，`scale.h = h*L/lambda_r`。

---

## File Structure

| 文件 | 操作 | 职责 |
|------|------|------|
| `example/热模块验证/jellyroll_vs_ring_thermal_compare.jl` | 修改 | 复制 Jellyroll scale 到 Ring；去掉热源桥接系数；可选增加验收标准与注释 |

---

## Chunk 1: 尺度对齐与热源修正

### Task 1: 将 Jellyroll 的 scale 复制到 Ring（在 SetCase 之前）

**Files:**
- Modify: `example/热模块验证/jellyroll_vs_ring_thermal_compare.jl`（`run_ring_simulation` 内，约 137–152 行）

- [ ] **Step 1: 在复制几何与热参数后、SetCase 前，复制 scale 中热相关字段**

在 `for prop in [:Rin, :Rout, :h, ...]` 循环之后、`opt = JuBat.Option()` 之前，增加：将 `param_dim`（Jellyroll）的 `scale` 中与热模型和归一化相关的字段复制到 `param_ring.scale`，使 Ring 与 Jellyroll 使用同一套参考尺度（从而 Bi、q_ref、L 一致）。

需复制的字段（与 `src/SetParams.jl` 中 ChooseCell 内设置及 `NormaliseParam` 使用一致）：`L`, `T_ref`, `t0`, `P_ref`, `lambda`, `q`, `h`。若 Ring 的 `param_dim.scale` 还参与其他量（如 `I_typ`, `phi`），为保持热一致性，可一并复制：`I_typ`, `phi`。

示例（在 `run_ring_simulation` 中，紧接几何复制循环后）：

```julia
# 与 Jellyroll 使用同一套无量纲尺度，使 Bi、q_ref、L 一致
for f in [:L, :T_ref, :t0, :P_ref, :lambda, :q, :h]
    setproperty!(param_ring.scale, f, getproperty(param_dim.scale, f))
end
# 可选：若 Ring 热求解或后续用到 I_typ/phi，也对齐
for f in [:I_typ, :phi]
    setproperty!(param_ring.scale, f, getproperty(param_dim.scale, f))
end
```

- [ ] **Step 2: 移除热源桥接系数，直接使用体热源**

将  
`q_bridge_scale = param_dim.cell.area`  
`q_mean_phys_ring = q_mean_phys .* q_bridge_scale`  
改为：  
`q_mean_phys_ring = q_mean_phys`  
并删除或注释掉对 `q_bridge_scale` 的打印。

- [ ] **Step 3: 运行验证脚本并检查曲线形态**

在项目根目录执行：

```bash
cd "d:\OneDrive\Desktop\Jubat For Cursor\JuBat"
julia --project=. "example/热模块验证/jellyroll_vs_ring_thermal_compare.jl"
```

检查：  
- Jellyroll 与 Ring 的温度曲线均应随时间的增加而单调上升（或近似单调），无“先上冲后明显下降”；  
- 两条曲线量级接近（温升在同一数量级）；  
- 若仍存在偏差，可再核对 `scale` 是否还有遗漏字段（如 `param_dim.scale` 的只读/派生字段）。

- [ ] **Step 4: 提交**

```bash
git add example/热模块验证/jellyroll_vs_ring_thermal_compare.jl
git commit -m "fix(thermal-verify): align Ring scale with Jellyroll and use direct volumetric heat source"
```

---

## Chunk 2: 文档与可维护性（可选）

### Task 2: 在脚本顶部或 run_ring_simulation 前增加简短说明

**Files:**
- Modify: `example/热模块验证/jellyroll_vs_ring_thermal_compare.jl`

- [ ] **Step 1: 添加注释说明尺度对齐与热源约定**

在 `run_ring_simulation` 中“创建 Ring 参数（几何与 Jellyroll 一致）”附近增加 1–2 行注释，说明：为与 Jellyroll 热模型可比，将 Jellyroll 的 scale（L, q, h 等）复制到 Ring，使无量纲化一致；热源为 Jellyroll 体积平均体热源 [W/m³]，直接用于 Ring 且用同一 q_ref 无量纲化。

- [ ] **Step 2: 提交（若做了 Chunk 2）**

```bash
git add example/热模块验证/jellyroll_vs_ring_thermal_compare.jl
git commit -m "docs(thermal-verify): comment scale alignment and heat source convention"
```

---

## Verification Checklist

- [ ] 脚本运行无报错。  
- [ ] Jellyroll 与 Ring 温度曲线均为单调上升（或近似），无圆环单独“上冲后下降”。  
- [ ] 两条温升量级接近（例如同属 0.x K 或 1.x K）。  
- [ ] 输出图与 CSV 与当前逻辑一致（仅数值与形态改善）。

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-03-14-jellyroll-ring-thermal-align.md`.  
Ready to execute with **superpowers:executing-plans** (batch of first 3 tasks, then report for feedback).
