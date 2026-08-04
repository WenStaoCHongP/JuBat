# CsvExport.jl Defensive Code Audit Plan

**Status:** ⬜ Pending | **Layer:** 5 后处理/导出 | **桶:** Mixed（Leave alone + category A 评估）

**Goal:** 7 处 try/catch 评估（CSV 写入失败容错）。CSV 写入失败应保留 @warn（用户磁盘/权限问题），但模式可统一为 helper。

## 现状（636 行）

| 区块 | 行号 | 用途 |
|---|---|---|
| `export_cycle_results_to_csv` | — | 主入口 |
| 7 处 try/catch | 98, 110, 127, 145, 156, 174, 190 | 每个 CSV 文件写入容错 |

**7 处 try/catch 模式**：
```julia
try
    write_csv(filepath, data)
catch e
    @warn "Failed to write <file>.csv" exception=e
end
```

**spec §7 评估**：**category B（资源失败容错）** — 不是 silent swallow，有 @warn + exception 记录。保留。

---

## Tasks

### Task 1: 评估 7 处 try/catch（保留 + 抽 helper）

**Files:** Modify `src/CsvExport.jl:98-200`

- [ ] **Step 1: 抽 helper 减重复**

在文件顶部添加：

```julia
function safe_write_csv(filepath::String, write_fn::Function, files_skipped::Vector{String})
    try
        write_fn()
    catch e
        @warn "Failed to write $(basename(filepath))" exception=e
        push!(files_skipped, basename(filepath))
    end
end
```

- [ ] **Step 2: 7 处调用改为 helper**

例如 line 98-103 改为：

```julia
safe_write_csv(joinpath(outdir, "cycle_summary.csv"),
    () -> write_csv_cycle_summary(joinpath(outdir, "cycle_summary.csv"), result),
    files_skipped)
```

- [ ] **Step 3: 跑主线 example**

Run: `julia example/czm_cycle_example.jl`
Expected: PASS（output/ 目录生成预期 CSV）

- [ ] **Step 4: Commit + baseline**

```bash
git add src/CsvExport.jl
git commit -m "refactor(csv): 7 处 try/catch 抽 safe_write_csv helper（spec §7 类目 B）"
echo "$(date +%F): CsvExport.jl 决策=[Consolidate: safe_write_csv helper]" >> Simplify/baseline.md
```

---

## Risk

低；helper 模式纯文本替换，行为等价。

## 不做的事

- 不删 try/catch（CSV 写入失败需要保留 @warn）
- 不动 CSV 文件格式
- 不合并多个 CSV 文件
