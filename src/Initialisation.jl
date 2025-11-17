function ModelInitialisation(case::Case)
    if isempty(case.opt.y0)
        if case.opt.model == "SPM"
            Nrn = case.mesh["negative particle"].nlen
            csn0 = ones(Float64, Nrn, 1) * case.param.NE.cs0
            Nrp = case.mesh["positive particle"].nlen
            csp0 = ones(Float64, Nrp, 1) * case.param.PE.cs0
            y0 = [csn0;  csp0]
        elseif case.opt.model == "SPMe"
            Nrn = case.mesh["negative particle"].nlen
            csn0 = ones(Float64, Nrn, 1) *  case.param.NE.cs0
            Nrp = case.mesh["positive particle"].nlen
            csp0 = ones(Float64, Nrp, 1) *  case.param.PE.cs0
            Ne = case.mesh["electrolyte"].nlen
            ce0 = ones(Float64, Ne, 1) *  case.param.EL.ce0
            y0 = [csn0;  csp0; ce0]
        elseif case.opt.model == "P2D"
            Nrn = case.mesh["negative particle"].nlen
            Nrp = case.mesh["positive particle"].nlen
            Ne = case.mesh["electrolyte"].nlen
            Nn = case.mesh["negative electrode"].nlen
            Np = case.mesh["positive electrode"].nlen
            csn0 = ones(Float64, Nrn, 1) * case.param.NE.cs0
            csp0 = ones(Float64, Nrp, 1) * case.param.PE.cs0
            ce0 = ones(Float64, Ne, 1) * case.param.EL.ce0
            phie0 = - ones(Float64, Ne, 1) * case.param.NE.U(case.param.NE.cs0)
            phis_p =  ones(Float64, Nn, 1) * case.param.PE.U(case.param.PE.cs0) .+ phie0[1] # guessed values are not used
            phis_n = zeros(Float64, Np, 1)
            y0 = [csn0;  csp0; ce0]
        else
            error( "Error: $(case.opt.model{1}) model has not been implemented!\n ")
        end
        if case.opt.thermalmodel == "lumped"
            y0 =[y0; case.param.cell.T0]
        elseif case.opt.thermalmodel == "distributed1D"
            # 预留：若未来加入 1D 分布式热网格，可在此追加 DOF
            # 当前无实现
        elseif case.opt.thermalmodel == "distributed2D"
            # 将二维分布式热的节点温度自由度追加到主状态向量，
            # 以便在时间推进中与电化学自由度一起求解（与 lumped 一致）。
            if haskey(case.mesh, "thermal2D")
                nT = case.mesh["thermal2D"].nlen
                T0_nodes = fill(case.param.cell.T0, nT)
                y0 = [y0; T0_nodes]
            end
        end
        if case.opt.model == "P2D"
            y0 =[y0; phis_n; phis_p; phie0]
        end
    else
        y0 = case.opt.y0 
    end
    return y0
end
