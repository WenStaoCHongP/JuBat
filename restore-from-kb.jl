"""
JuBat ↔ lib-knowledge-base 归档恢复工具
（规格见已归档规格 raw/repos/jubat/docs-archive：2026-09-09 JuBat 瘦身与知识库收录）

用法（仓库根执行）：
  julia restore-from-kb.jl --what test|tools|md|all [--kb <路径>]   # KB → 工作区
  julia restore-from-kb.jl --sync --what test|tools|all             # 工作区 → KB（覆盖同路径）

--what md 在 KB 侧对应 raw/repos/jubat/md-raw/。覆盖前打印将被覆盖的文件数。
"""
const DEFAULT_KB = raw"D:\OneDrive\Desktop\lib知识库\lib-knowledge-base"

kb_name(what) = what == "md" ? "md-raw" : what

function parse_args(args)
    what = ""; kb = DEFAULT_KB; sync = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--what" && i < length(args)
            i += 1; what = args[i]
        elseif a == "--kb" && i < length(args)
            i += 1; kb = args[i]
        elseif a == "--sync"
            sync = true
        else
            error("未知参数：$a（用法见文件头注释）")
        end
        i += 1
    end
    what in ("test", "tools", "md", "all") ||
        error("--what 必须是 test|tools|md|all，got $what")
    sync && what == "md" && error("--sync 不支持 md（wiki 为编译产物，单向）")
    sync && what == "all" && error("--sync 请逐项指定 test/tools")
    return what, kb, sync
end

countfiles(p) = isdir(p) ? count(true for (_, _, files) in walkdir(p) for _ in files) : 0

function transfer(src, dst)
    isdir(src) || error("源目录不存在：$src")
    if isdir(dst)
        println("将覆盖已存在的目标：$dst（现有 $(countfiles(dst)) 文件 ← 源 $(countfiles(src)) 文件）")
    end
    cp(src, dst; force=true)
    println("完成：$src → $dst（$(countfiles(dst)) 文件）")
end

function main()
    what, kb, sync = parse_args(ARGS)
    rawroot = joinpath(kb, "raw", "repos", "jubat")
    items = what == "all" ? ["test", "tools", "md"] : [what]
    for w in items
        if sync
            transfer(joinpath(@__DIR__, w), joinpath(rawroot, kb_name(w)))
        else
            transfer(joinpath(rawroot, kb_name(w)), joinpath(@__DIR__, w))
        end
    end
end
main()
