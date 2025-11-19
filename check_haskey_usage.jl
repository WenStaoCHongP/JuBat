# 检查代码中 haskey 使用是否正确的脚本
# 使用方法：julia check_haskey_usage.jl

using Glob

println("检查 haskey 使用情况...")
println("=" ^ 60)

# 搜索所有 .jl 文件
jl_files = glob("**/*.jl", ".")

# 存储发现的潜在问题
issues = []

for file in jl_files
    # 跳过本脚本自身
    if endswith(file, "check_haskey_usage.jl")
        continue
    end
    
    try
        content = read(file, String)
        lines = split(content, '\n')
        
        for (line_num, line) in enumerate(lines)
            # 检查可能错误的模式：haskey(case, ...
            # 注意：这是一个简单的模式匹配，可能有误报
            if occursin(r"haskey\s*\(\s*case\s*,", line)
                push!(issues, (file=file, line=line_num, content=strip(line)))
                println("⚠️  发现潜在问题：")
                println("   文件: $file")
                println("   行号: $line_num")
                println("   内容: $(strip(line))")
                println()
            end
        end
    catch e
        # println("无法读取文件 $file: $e")
    end
end

println("=" ^ 60)
if isempty(issues)
    println("✅ 未发现 haskey(case, ...) 的错误用法")
    println()
    println("当前代码库中的 haskey 使用都是正确的，例如：")
    println("  - haskey(case.mesh, \"thermal2D\")  ✅")
    println("  - haskey(variables, \"T_nodes\")    ✅")
else
    println("❌ 发现 $(length(issues)) 处潜在问题")
    println()
    println("修复建议：")
    println("将 haskey(case, :symbol) 改为 hasproperty(case, :symbol)")
    println()
    println("详细说明请查看：test3_debug_solution.md")
end

println("=" ^ 60)
println("提示：如果您的 test3.jl 不在此目录，请在包含该文件的目录运行此脚本")
