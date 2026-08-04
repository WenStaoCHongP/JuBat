# parameters/Jellyroll.jl Audit Plan

**Status:** ⬜ Pending
**Layer:** 1 - 参数与变量
**桶:** Leave alone（主线必保）

**Goal:** 审查 `parameters/Jellyroll.jl`；确认 TODO 占位在 md 文档跟踪，不在代码中清理。

---

## Audit Checklist

- [ ] **Step 1: 通读 `Jellyroll.jl` (202 行)**，列出所有 `TODO 用户提供实测值`（spec §7.4 类别 D）

- [ ] **Step 2: 把 TODO 清单复制到 `md/01_参数定义与归一化.md`**（参数文档跟踪）

- [ ] **Step 3: 跑 `ChooseCell("Jellyroll")` 验证**

Run: `julia -e 'include("src/JuBat.jl"); using .JuBat; ChooseCell("Jellyroll"); println("OK")'`
Expected: 无错（即使 `@warn` 提示 E_coat/cohesive 缺失也可接受，非本次范围）

- [ ] **Step 4: baseline 记录**

```bash
echo "$(date +%F): parameters/Jellyroll.jl 已审查，TODO 转交 md/01" >> Simplify/baseline.md
```

---

## Result

无代码修改；TODO 数据缺口转交参数文档跟踪。
