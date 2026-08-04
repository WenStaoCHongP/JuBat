# VariableKeys.jl Implementation Plan

**Status:** ⬜ Pending
**Layer:** 1 - 参数与变量
**桶:** 新建文件（D2 重复簇根因修复）
**依赖:** 无（其他文件依赖此文件）
**阻塞:** Variables.md、CallModel.md、CouplingState.md 必须等本文件完成

**Goal:** 新建 `src/VariableKeys.jl`，集中所有 `variables["..."]` 字符串键名为常量，消除 Variables.jl / CallModel.jl / CouplingState.jl 三处硬编码键表（spec D2）。

**Architecture:** 一个 `const` 字典 + 多个 `const` 元组，按"全局变量键 / 单元变量键 / 单元结果键"分类。其他文件引用常量而非字面量。

**Tech Stack:** Julia 1.x；无外部依赖。

---

## Files

- Create: `src/VariableKeys.jl`
- Modify: `src/JuBat.jl:4`（include 顺序最早）

---

## Tasks

### Task 1: 抽取键名清单

**Files:**
- Read: `src/Variables.jl:1-156`（`StandardVariables` 全部键）
- Read: `src/Variables.jl:157-235`（`create_element_workspace` 全部键）
- Read: `src/CallModel.jl:247-268`（`copy_element_results` 全部键）

- [ ] **Step 1.1: 列举 StandardVariables 键**

读 `Variables.jl:1-156`，把所有 `"..."=>...` 的键名抽取到列表 `GLOBAL_KEYS`。

- [ ] **Step 1.2: 列举 create_element_workspace 键**

读 `Variables.jl:157-235`，抽取到 `ELEMENT_WORKSPACE_KEYS`。

- [ ] **Step 1.3: 列举 copy_element_results 键**

读 `CallModel.jl:247-268`，抽取到 `ELEMENT_RESULT_KEYS`。

- [ ] **Step 1.4: 标注重叠键**

三个清单中重叠的键（如 `"temperature"`、`"cell voltage"`）标出，验证语义一致后合并到 `COMMON_KEYS`。

- [ ] **Step 1.5: Commit 键清单（暂存到草稿）**

```bash
git add -A
git commit -m "chore(VariableKeys): draft key inventory"
```

### Task 2: 写入 `src/VariableKeys.jl`

**Files:**
- Create: `src/VariableKeys.jl`

- [ ] **Step 2.1: 写文件骨架**

```julia
# src/VariableKeys.jl
# 集中所有 variables["..."] 字符串键常量
# 详见 docs/superpowers/specs/2026-08-04-src-simplification-design.md D2

module VariableKeys

# 全局变量键（StandardVariables 用）
const GLOBAL_KEYS = (
    # === 状态提取 ===
    "negative particle lithium concentration",
    "positive particle lithium concentration",
    "negative particle surface lithium concentration",
    "positive particle surface lithium concentration",
    # === 电压/电流 ===
    "cell voltage",
    "negative electrode overpotential",
    # ... Task 1 抽取的全部键
)

# 单元 workspace 键（create_element_workspace 用）
const ELEMENT_WORKSPACE_KEYS = (
    "negative electrode overpotential",
    "positive electrode overpotential",
    # ... Task 1 抽取
)

# 单元结果键（copy_element_results 用）
const ELEMENT_RESULT_KEYS = (
    "negative electrode overpotential",
    "positive electrode overpotential",
    "cell voltage",
    # ... Task 1 抽取
)

end  # module
```

- [ ] **Step 2.2: 完整填充所有键**

按 Task 1 的清单逐字面量填入，**保留原拼写与空格**（任何键名拼写差异会破坏 477 处 `variables[...]` 访问）。

- [ ] **Step 2.3: 在 `JuBat.jl:4` 之前 include**

修改 `src/JuBat.jl`，在第一个 `include("Option.jl")` 之前加：

```julia
include("VariableKeys.jl")  # 键常量定义（D2 根因修复）
```

- [ ] **Step 2.4: 验证模块加载**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; println(VariableKeys.GLOBAL_KEYS[1])'`
Expected: 打印 `"negative particle lithium concentration"`

- [ ] **Step 2.5: Commit**

```bash
git add src/VariableKeys.jl src/JuBat.jl
git commit -m "feat(VariableKeys): 新建集中化变量键常量模块（D2 根因修复）"
```

---

## Validation

- [ ] `using .JuBat` 无报错
- [ ] `VariableKeys.GLOBAL_KEYS` 长度等于原 `StandardVariables` 键数
- [ ] `VariableKeys.ELEMENT_RESULT_KEYS` 长度等于 `copy_element_results` 键数（应为 16）

## Risk

- **键名拼写偏差**：任何拼写差异会破坏 17 文件 477 处 `variables["..."]` 访问。**缓解**：Step 2.4 验证；后续 Variables.md plan 加 characterization test
- **模块加载顺序**：必须在所有使用键的文件之前 include。已在 `JuBat.jl:4` 之前位置保证
