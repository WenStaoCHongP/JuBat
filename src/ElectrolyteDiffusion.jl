function  ElectrolyteDiffusion(param::Params, mesh::Mesh, mlen::Int64, variables::Dict{String, Union{Array{Float64},Float64}})
    """
        Generate the system equatiuon of M du/dt=Ku+F
        inputs = electrode -- eledtrode parameter
                 mesh -- mesh
                 mlen -- length of unknowns
                 variables -- dictionary of variables
        outputs = M -- mass matrix
                  K -- stiffness matrix
    """
    ce_n_gs = variables["electrolyte lithium concentration at negative electrode Gauss point"] 
    ce_p_gs = variables["electrolyte lithium concentration at positive electrode Gauss point"] 
    ce_sp_gs = variables["electrolyte lithium concentration at separator Gauss point"]
    M = electrolyte_diffusion_mass(param, mesh, mlen, size(ce_n_gs), size(ce_sp_gs), size(ce_p_gs))
    K = electrolyte_diffusion_stiffness(param, mesh, mlen, variables)
    return M, K
end

function electrolyte_diffusion_mass(param::Params, mesh::Mesh, mlen::Int64, shape_ne, shape_sp, shape_pe)
    Vi = mesh.element[mesh.gs.ele,:]
    coeff_ne = ones(shape_ne) * param.NE.eps
    coeff_sp = ones(shape_sp) * param.SP.eps
    coeff_pe = ones(shape_pe) * param.PE.eps
    coeff = [coeff_ne; coeff_sp; coeff_pe] .* mesh.gs.weight .* mesh.gs.detJ

    return Assemble(Vi, Vi, mesh.gs.Ni, mesh.gs.Ni, coeff, mlen)
end

function electrolyte_diffusion_stiffness(param::Params, mesh::Mesh, mlen::Int64, variables::Dict{String, Union{Array{Float64},Float64}})
    Vi = mesh.element[mesh.gs.ele,:]
    ce_n_gs = variables["electrolyte lithium concentration at negative electrode Gauss point"]
    ce_p_gs = variables["electrolyte lithium concentration at positive electrode Gauss point"]
    ce_sp_gs = variables["electrolyte lithium concentration at separator Gauss point"]
    T = variables["temperature"]
    De_ne =  param.EL.De(ce_n_gs,T) * param.NE.eps ^ param.NE.brugg
    De_sp =  param.EL.De(ce_sp_gs,T) * param.SP.eps ^ param.SP.brugg
    De_pe =  param.EL.De(ce_p_gs,T) * param.PE.eps ^ param.PE.brugg
    De_eff =  [De_ne; De_sp; De_pe] * Arrhenius(param.EL.Eac_D, T)
    coeff = - De_eff .* mesh.gs.weight .* mesh.gs.detJ
    return Assemble(Vi, Vi, mesh.gs.dNidx, mesh.gs.dNidx, coeff, mlen)
end
