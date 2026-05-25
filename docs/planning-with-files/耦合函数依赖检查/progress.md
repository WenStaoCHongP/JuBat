# 耦合函数依赖检查 - 进度

## 状态: 已完成

- [x] 读取所有核心 src 文件
- [x] 梳理函数调用关系
- [x] 撰写 findings.md
  - [x] Part A: 执行流程总览 (A1-A5)
    - [x] A1 初始化阶段
    - [x] A2 求解阶段主循环
    - [x] A3 CallModel_MultiSPMe 内部流程
    - [x] A4 CZM 损伤更新流程
    - [x] A5 数据流概要
  - [x] Part B: 模块函数清单 (B1-B14)
    - [x] B1 SetParams.jl
    - [x] B2 Option.jl
    - [x] B3 SetCase.jl
    - [x] B4 CouplingState.jl
    - [x] B5 CallModel.jl
    - [x] B6 SPMe.jl
    - [x] B7 Mechanical.jl
    - [x] B8 Parallelsolution.jl
    - [x] B9 ThermalDistributed.jl
    - [x] B10 czm.jl
    - [x] B11 CzmSolve.jl
    - [x] B12 Solve.jl
    - [x] B13 Jellyrollmodel.jl
    - [x] B14 CycleSolver.jl
  - [x] 附录: 完整开关→函数映射表
