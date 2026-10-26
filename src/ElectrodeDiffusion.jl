function  ElectrodeDiffusion(electrode::Electrode, mesh::Mesh, mlen::Int64, c::Matrix{Float64}=zeros(Float64), theta::Float64=0.0, T::Float64=1.0)
    """
        Generate the system equatiuon of M du/dt=Ku+F
        inputs = electrode -- eledtrode parameter
                 mesh -- mesh
                 mlen -- length of unknowns
                 T -- temperature
        outputs = M -- mass matrix
                  K -- stiffness matrix
    """
    M = electrode_diffusion_mass(mesh, mlen)
    K = electrode_diffusion_stiffness(electrode, mesh, mlen, c, theta, T)
    return M, K
end

function electrode_diffusion_mass(mesh::Mesh, mlen::Int64)
    Vi = mesh.element[mesh.gs.ele,:]
    coeff = mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ
    return Assemble(Vi, Vi, mesh.gs.Ni, mesh.gs.Ni, coeff, mlen)
end

function electrode_diffusion_stiffness(electrode::Electrode, mesh::Mesh, mlen::Int64, c::Matrix{Float64}, theta::Float64, T::Float64=1.0)
    Ds_eff =  electrode.Ds * Arrhenius(electrode.Eac_D, T).* (1 .+ (theta.*c))
    return electrode_diffusion_stiffness(mesh, mlen, Ds_eff)
end

function electrode_diffusion_stiffness(mesh::Mesh, mlen::Int64, Ds_eff)
    Vi = mesh.element[mesh.gs.ele,:]
    coeff = - Ds_eff .* mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ
    return Assemble(Vi, Vi, mesh.gs.dNidx, mesh.gs.dNidx, coeff, mlen)
end
