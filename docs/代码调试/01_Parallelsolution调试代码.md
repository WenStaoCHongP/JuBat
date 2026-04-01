# Parallelsolution.jl 调试代码记录

> 文件: `src/Parallelsolution.jl`
> 状态: A (新增文件)
> 移除行数: ~100 行

---

## 1. 调试函数 (完整删除)

### 1.1 debug_check_prefactors (原 lines 5-30)

**功能**: 检查电化学预因子是否包含 NaN/Inf，打印详细诊断信息。

```julia
# 调试输出函数
"""检查并报告 NaN/Inf 值（调试用）"""
function debug_check_prefactors(prefactor_n, prefactor_p, csn_av, csp_av, u_n_ref_val, u_p_ref_val, du_n_dT_val, du_p_dT_val, c_sigma, cn_surf, cp_surf, ce_n_gs, ce_p_gs)
	has_nan = !isfinite(prefactor_n) || !isfinite(prefactor_p) || !isfinite(csn_av) || !isfinite(csp_av) || !isfinite(u_n_ref_val) || !isfinite(u_p_ref_val) || !isfinite(c_sigma)

	if !has_nan
		return false
	end

	println("\n" * "="^80)
	println("❌ [DEBUG] 预计算值包含 NaN/Inf")
	println("="^80)
	println("📊 预因子: prefactor_n=$prefactor_n, prefactor_p=$prefactor_p")
	println("📊 平均浓度: csn_av=$csn_av, csp_av=$csp_av")
	println("📊 开路电位: u_n=$u_n_ref_val, u_p=$u_p_ref_val")
	println("📊 温度系数: du_n_dT=$du_n_dT_val, du_p_dT=$du_p_dT_val")
	println("📊 电导: c_sigma=$c_sigma")
	println("\n💡 输入浓度前5个值:")
	println("  cn_surf: $(cn_surf[1:min(5,end)])")
	println("  cp_surf: $(cp_surf[1:min(5,end)])")
	println("  ce_n_gs: $(ce_n_gs[1:min(5,end)])")
	println("  ce_p_gs: $(ce_p_gs[1:min(5,end)])")
	println("="^80 * "\n")

	return true
end
```

**恢复**: 粘贴到 `compute_prefactors()` 函数之后。

---

### 1.2 debug_check_coefficients (原 lines 32-52)

**功能**: 检查单元电化学系数是否异常，打印诊断信息。

```julia
"""检查单元系数是否有效"""
function debug_check_coefficients(e, has_nan_prefactor, C1, C2, alpha_p, alpha_n, C5, T_e, j0_n, j0_p, u_n_val_T, u_p_val_T)
	if has_nan_prefactor || e > 1
		return  # 只打印第一个异常单元，且预计算值正常时
	end

	has_nan_coeff = !isfinite(C1) || !isfinite(C2) || !isfinite(alpha_p) || !isfinite(alpha_n) || !isfinite(C5)

	if !has_nan_coeff
		return
	end

	println("\n" * "="^80)
	println("❌ [DEBUG] 单元 $e 系数异常")
	println("="^80)
	println("📊 T_e=$T_e, j0_n=$j0_n, j0_p=$j0_p")
	println("📊 u_n(T)=$u_n_val_T, u_p(T)=$u_p_val_T")
	println("📊 系数: C1=$C1, C2=$C2")
	println("   alpha_p=$alpha_p, alpha_n=$alpha_n, C5=$C5")
	println("="^80 * "\n")
end
```

**恢复**: 粘贴到 `debug_check_prefactors()` 之后。

---

### 1.3 debug_check_initial_voltage (原 lines 54-76)

**功能**: 检查初始电压是否有效，打印异常分支诊断。

```julia
"""检查初始电压是否有效"""
function debug_check_initial_voltage(has_nan_prefactor, V, V_branches, I_e, coeffs, I_total, ne)
	if has_nan_prefactor || isfinite(V)
		return
	end

	println("\n" * "="^80)
	println("❌ [DEBUG] 初始电压异常: V=$V")
	println("="^80)
	println("  I_total=$I_total, ne=$ne")
	println("  前3个异常单元:")

	count = 0
	for e in 1:ne
		if !isfinite(V_branches[e]) && count < 3
			count += 1
			println("  单元 $e: V=$(V_branches[e]), I=$(I_e[e])")
			println("    C1=$(coeffs[e].C1), C2=$(coeffs[e].C2)")
			println("    α_p=$(coeffs[e].alpha_p), α_n=$(coeffs[e].alpha_n)")
		end
	end
	println("="^80 * "\n")
end
```

**恢复**: 粘贴到 `debug_check_coefficients()` 之后。

---

## 2. 调试模式参数 (移除参数和调用)

### 2.1 compute_element_coefficients 的 debug_mode 参数 (原 line 119)

**原签名**:
```julia
function compute_element_coefficients(e, T_e, param, prefactors, T_ref, debug_mode=false)
```

**改为**:
```julia
function compute_element_coefficients(e, T_e, param, prefactors, T_ref)
```

### 2.2 compute_element_coefficients 内的 debug_mode 条件块 (原 lines 143-146)

**原代码**:
```julia
	# 调试检查（仅第一个异常单元）
	if debug_mode
		debug_check_coefficients(e, false, C1, C2, alpha_p, alpha_n, C5, T_e, j0_n, j0_p, u_n, u_p)
	end
```

**整块删除**。

### 2.3 compute_all_coefficients 的 debug_mode 参数 (原 line 152)

**原签名**:
```julia
function compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref, debug_mode=false)
```

**改为**:
```julia
function compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref)
```

### 2.4 compute_all_coefficients 内传递 debug_mode (原 line 155)

**原代码**:
```julia
		coeffs[e] = compute_element_coefficients(e, Te_prev[e], param, prefactors, T_ref, debug_mode)
```

**改为**:
```julia
		coeffs[e] = compute_element_coefficients(e, Te_prev[e], param, prefactors, T_ref)
```

---

## 3. solve_branch_currents_newton 中的调试调用 (原 lines 478-491, 534-536)

### 3.1 debug_mode 变量和 prefactor 检查 (原 lines 478-487)

**原代码**:
```julia
	# 调试：检查预因子
	debug_mode = case.opt.debug_coupling
	has_nan_prefactor = debug_check_prefactors(
		prefactors.prefactor_n, prefactors.prefactor_p,
		prefactors.csn_av, prefactors.csp_av,
		prefactors.u_n_ref_val, prefactors.u_p_ref_val,
		prefactors.du_n_dT_val, prefactors.du_p_dT_val,
		prefactors.c_sigma,
		prefactors.cn_surf, prefactors.cp_surf,
		prefactors.ce_n_gs, prefactors.ce_p_gs
	)
```

**整块替换为**:
```julia
	# 预因子已在 compute_prefactors 中计算
```

### 3.2 compute_all_coefficients 调用中的 debug_mode (原 line 491)

**原代码**:
```julia
	coeffs = compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref, debug_mode)
```

**改为**:
```julia
	coeffs = compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref)
```

### 3.3 V_branches_all 和 debug_check_initial_voltage 调用 (原 lines 534-536)

**原代码**:
```julia
	# 调试：检查初始电压
	V_branches_all = [branch_voltage(coeffs[e], I_e[e]) for e in 1:ne]
	debug_check_initial_voltage(has_nan_prefactor, V, V_branches_all, I_e, coeffs, I_total, ne)
```

**整块删除** (V_branches_all 与上方 V_branches 计算重复)。

---

## 4. 恢复依赖

恢复调试代码时需同步恢复:
1. `src/Option.jl` 中的 `debug_coupling::Bool = false` 字段（保留未删除）
2. 以上 3 个调试函数粘贴回原位置
3. `compute_element_coefficients` 和 `compute_all_coefficients` 添加 `debug_mode` 参数
4. `solve_branch_currents_newton` 中恢复调试调用代码块
