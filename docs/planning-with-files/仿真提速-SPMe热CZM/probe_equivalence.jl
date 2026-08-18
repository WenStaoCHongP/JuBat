# 等价性 1：稀疏方阵 lu(A)\b 与 A\b 是否逐位一致（决定 Task 6 分解因子缓存可行性）
# 等价性 2：同 pattern 稀疏矩阵手写同序加法与 SparseArrays + 是否逐位一致（决定 Task 7 D1 可行性）
using SparseArrays, LinearAlgebra, Random
Random.seed!(20260819)

ok1 = true
for trial in 1:20
    n = rand(200:800)
    A = sprand(n, n, 0.005); A = A + A' + 2.0I   # 对称正定近似，模拟 M - θKdt
    b = randn(n)
    if (lu(A) \ b) != (A \ b); ok1 = false; end
end
println("lu(A)\\b == A\\b : ", ok1)

ok2 = true
for trial in 1:20
    n = rand(400:1200)
    Kb = sprand(n, n, 0.004)                     # 模拟 K_bulk
    Kc = sprand(n, n, 0.0008)                    # 模拟 K_coh（更稀疏）
    Kt = Kb + Kc                                 # 并集 pattern，只算一次
    rv = rowvals(Kt)
    manual = copy(nonzeros(Kt))
    for col in 1:size(Kt, 2)                     # 对 Kt 每个存储位置恰一次 a+b 加法
        for k in Kt.colptr[col]:(Kt.colptr[col+1]-1)
            i = rv[k]
            manual[k] = Kb[i, col] + Kc[i, col]
        end
    end
    Kman = SparseMatrixCSC(size(Kt,1), size(Kt,2), copy(Kt.colptr), copy(rv), manual)
    if Kman != Kt; global ok2 = false; end
end
println("手写同序加法 == SparseArrays + : ", ok2)
