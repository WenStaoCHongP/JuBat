# parameters/Enertech.jl Dead Code Verification Plan

**Status:** ✅ Completed（保留）
**Layer:** 1 - 参数与变量
**桶:** Delete 候选

**Goal:** grep 验证 `ChooseCell("Enertech")` 是否有调用者；无则删除整个文件 + 移除 JuBat.jl 中的相关 export。

---

## Tasks

### Task 1: grep 验证

- [ ] **Step 1.1: 全仓库搜调用点**

```bash
grep -rn 'ChooseCell("Enertech"' --include="*.jl" .
```

记录匹配数与位置。

- [ ] **Step 1.2: 搜其他可能引用**

```bash
grep -rn "Enertech" --include="*.jl" --include="*.md" --include="*.toml" .
```

包括 import、documentation、Project.toml 等。

### Task 2: 决策与执行

- [ ] **Step 2.1: 若 Step 1.1 仅匹配 `SetParams.jl:ChooseCell` 定义本身**

→ 整文件删除：

```bash
git rm src/parameters/Enertech.jl
```

- [ ] **Step 2.2: 若 SetParams.jl 的 `ChooseCell` 有 `"Enertech"` 分支**

修改 `src/SetParams.jl:249-349`，移除 `"Enertech" => ...` 分支。

- [ ] **Step 2.3: 移除 JuBat.jl 中相关 export（如有）**

```bash
grep -n "Enertech" src/JuBat.jl
```

如有，删除对应行。

- [ ] **Step 2.4: 跑主线 example 验证未破坏**

Run: `julia example/minimal_example.jl`
Expected: PASS（Jellyroll 主线不依赖 Enertech）

- [ ] **Step 2.5: Commit**

```bash
git add -A
git commit -m "chore(parameters): 删除未使用的 Enertech 参数集（dead code）"
```

- [ ] **Step 2.6: 若 Step 1.1 有任何外部调用者**

→ 改桶为 Leave alone；在 baseline.md 标"保留——被 X 使用"。

### Task 3: baseline 更新

- [ ] **Step 3.1: 记录**

```bash
echo "$(date +%F): parameters/Enertech.jl 决策=[删除/保留]" >> Simplify/baseline.md
```

---

## Risk

- 低；参数集删除影响范围仅限 `ChooseCell` 调用者
- **必须先 grep 再删**，不能凭直觉

## Execution Result (2026-08-05)

- 仓库示例未直接调用，但 `ChooseCell` 文档明确将 `Enertech` 列为支持参数集，属于公开字符串 API。
- 仅凭仓库内 grep 不足以证明外部用户不存在；按评审后的公共 API 规则保留。
