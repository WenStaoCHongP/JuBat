using SparseArrays: SparseMatrixCSC

mutable struct GaussPoint
    x::Array{Float64}   # physical coordinates
    xi::Array{Float64}  # local coordinates
    weight::Vector{Float64}
    detJ::Vector{Float64}
    ele::Vector{Int64}  # element information
    Ni::Array{Float64}  # shape functions
    dNidx::Array{Float64}   # derivatives of shape functions
    order::Int64    # the order of Gauss quadrature
end

mutable struct Mesh
    type::String
    dimension::Int64
    node::Array{Float64}
    nlen::Int64 # the length of nodes
    element::Array{Int64}
    gs::GaussPoint # Gauss points
end

# Abstract CZM types are used by CohesiveMesh and defined here to avoid
# include-order issues with czm.jl.
abstract type AbstractCohesiveElement end
abstract type AbstractDamageState end

"""
    CzmSubmesh

径向 8 层的 CZM 机械子网格。周向节点继承粗热网格，
通过 thermal_elem_map 与 thermal_to_czm 矩阵耦合。

# 字段
- `mesh`: 细化 Q4 网格
- `material_type`: 每个单元的材料类型（:PE / :PCC / :SP / :NE / :NCC）
- `winding_turn`: 卷绕圈号（从内到外 1, 2, ...）
- `thermal_elem_map`: 每个 CZM 单元 → 对应的粗热单元 id
- `phi_pairs`: 力学外螺旋 θ 节点 → 下一匝内螺旋 θ+2π 节点对

定义在 SetMesh.jl（而非 czm.jl）以避免 include 顺序问题：
CohesiveMesh 引用 CzmSubmesh，而 SetMesh.jl 在 czm.jl 之前被 include。
"""
struct CzmSubmesh
    mesh::Mesh                              # 细化 Q4 网格
    material_type::Vector{Symbol}           # :PE / :PCC / :SP / :NE / :NCC
    winding_turn::Vector{Int}               # 卷绕圈号（从内到外 1, 2, ...）
    thermal_elem_map::Vector{Int}           # 每个 CZM 单元 → 对应的粗热单元 id
    phi_pairs::Vector{Tuple{Int,Int}}        # (outer_node_at_θ, inner_node_at_θ+2π)
end

mutable struct CohesiveMesh
    bulk_mesh::Mesh                           # 原始固体网格
    node::Matrix{Float64}                     # 扩展后节点坐标
    nnode::Int64                              # 总节点数
    bulk_element::Matrix{Int64}               # 更新后的固体单元连接关系
    cohesive_elements::Vector{AbstractCohesiveElement} # 内聚力单元
    n_cohesive::Int64                         # 内聚力单元数
    n_layers::Int64                           # 遗留字段名：实际保存 cohesive 本构/材料类型数（2），不是物理层数、真实面数或单元数
    node_map::Dict{Int64, Vector{Int64}}      # 原节点 → [分层后的节点们]
    interface_nodes::Vector{Vector{Tuple{Int64,Int64}}} # 每个界面的节点对
    damage_states::Vector{AbstractDamageState} # 损伤状态
    czm_submesh::Union{Nothing, CzmSubmesh}                       # v5 新增：细化 CZM 子网格
    thermal_to_czm::Union{Nothing, SparseMatrixCSC{Float64, Int}} # v5 新增：粗热 → 细 CZM 插值矩阵
    cohesive_to_thermal::Union{Nothing, Vector{Int}}              # v5 新增：CZM 单元 → 粗热单元 id 反向映射

    # 内部构造函数（空初始化）
    function CohesiveMesh()
        new(Mesh("Q4", 2, zeros(0,2), 0, zeros(Int64,0,4),
            GaussPoint(zeros(0,2), zeros(0,2), zeros(0), zeros(0), zeros(Int64,0), zeros(0,4), zeros(0,8), 2)),
            zeros(0, 2), 0, zeros(Int64, 0, 4),
            AbstractCohesiveElement[], 0, 0, Dict{Int64, Vector{Int64}}(),
            Vector{Vector{Tuple{Int64,Int64}}}(), AbstractDamageState[],
            nothing, nothing, nothing)  # 新字段默认 nothing
    end
end

function SetMesh(domain::Any, num::Any, type::String, gsorder::Int64=4)
"""
    A function to set up mesh
     inputs are 'domain, num, type, gsorder'
     'domain' is an array for the problem domain information
     'num' is a vector for the number of elements
     'type' is a string for the element type
     'gsorder' is int for the order of Gauss quadrature, with default value 4
     e.g. mesh =  SetMesh([0,1,2],[5,5],'L2')
     There should be more sophisticated methods to build a mesh, to be finished
"""

    if type in ["L2", "L3"]
        mesh = Mesh1D(domain, num, type, gsorder)
    elseif type == "Q4"
        mesh = Mesh2D(domain, num, type, gsorder)
    elseif type == "COH2D4"
        mesh = Mesh2D(domain, num, type, gsorder)
    else
        error("Error: element type $type has not been implemented!\n")
    end
    return mesh
end

function Mesh2D(domain::Vector{Float64}, num::Any, type::String="Q4", gsorder::Int64=2)
    """
        构建规则矩形区域的二维Q4网格
        domain = [x0, x1, y0, y1]
        num = [nx, ny] 或 Tuple(nx, ny)
    """
    dim = 2
    @assert length(domain) == 4 "domain for Q4 must be [x0,x1,y0,y1]"
    if isa(num, Tuple) || isa(num, NTuple{2,Int})
        nx, ny = num
    else
        @assert length(num) == 2 "num for Q4 must have length 2"
        nx = Int(num[1]); ny = Int(num[2])
    end
    x0, x1, y0, y1 = domain
    dx = (x1 - x0)/nx
    dy = (y1 - y0)/ny
    # nodes: (nx+1)*(ny+1)
    nnx = nx + 1; nny = ny + 1
    nnode = nnx * nny
    node = zeros(Float64, nnode, dim)
    # ordering: loop j (y) outer, i (x) inner; index(i,j) = j*nnx + i + 1 ; i,j start at 0
    for j in 0:ny
        for i in 0:nx
            idx = j*nnx + i + 1
            node[idx,1] = x0 + i*dx
            node[idx,2] = y0 + j*dy
        end
    end
    # elements: nx * ny, node ordering consistent with Q4 reference (1:(i,j),2:(i+1,j),3:(i+1,j+1),4:(i,j+1))
    ne = nx * ny
    element = zeros(Int64, ne, 4)
    e = 0
    for j in 0:ny-1
        for i in 0:nx-1
            e += 1
            n1 = j*nnx + i + 1
            n2 = n1 + 1
            n3 = n2 + nnx
            n4 = n1 + nnx
            # 逐元素赋值避免 Julia 对行切片的单值广播限制
            element[e,1] = n1
            element[e,2] = n2
            element[e,3] = n3
            element[e,4] = n4
        end
    end
    gs = GetGS(element, node, gsorder, type)
    mesh = Mesh(type, dim, node, nnode, element, gs)
    return mesh
end

function Mesh1D(domain::Vector{Float64}, num::Any, type::String= "L2", gsorder::Int64= 4)
# to build 1D mesh and its gauss points
    dim = 1
    element_number = round(Int, sum(num))
    if type == "L2"
        ele_node = 2
    elseif type == "L3"
        ele_node = 3
    end
    node = zeros(Float64, element_number * (ele_node-1) + 1, dim)
    v = 0
    for i = 1:length(domain) - 1
        dx = 1/num[i]/(ele_node - 1)
        temp = collect(0:dx:1) .* (domain[i+1] - domain[i]) .+ domain[i]
        len = num[i] * (ele_node - 1)
        if i == 1
            node[v + 1:v + len + 1, 1] = temp
            v = v + len + 1
        else
            node[v + 1:v + len, 1] = temp[2:end,1]
            v = v + len
        end
    end
    nlen = deepcopy(v)
    element = zeros(Int64, element_number, ele_node)
    v = 0 
    ele = 0
    for i in num
        for j = 1:i
            ele = ele + 1
            element[ele, 1:ele_node] = v + 1:v + ele_node
            v = v + ele_node - 1
        end
    end
    gs = GetGS(element, node, gsorder, type)
    mesh = Mesh(type, dim, node, nlen, element, gs)
    return mesh
end

function PickElement(mesh::Mesh, v::Vector{Int64})
    """
        pick elements from a mesh
        Inputs = mesh::Mesh
                v::Vector{Inveger}, the index of elements to be picked up
        outputs = mesh1::Mesh
            mesh1 = mesh(v), including the collection of Gaussian points
        the new mesh will resort the index of node, element and Gauss points 
    """
    type = mesh.type
    dim = mesh.dimension
    element = deepcopy(mesh.element[v,:])
    node_pool = sort(unique(reshape(element,:,1)))
    node_pool_pair = zeros(Int64, maximum(node_pool))
    nlen = length(node_pool)
    node_pool_pair[node_pool] = collect(1:nlen)
    node = mesh.node[node_pool,:]
    element = node_pool_pair[element]
    gsorder = mesh.gs.order
    gs = deepcopy(mesh.gs)
    len = gsorder ^ dim
    v_gs = zeros(Int64, length(v) * len)
    gs.ele = zeros(Int64, length(v) * len)
    for i = 1:len
        v_gs[i:len:length(v) * len] = (v .- 1) .* len .+ i
        gs.ele[i:len:length(v) * len] = 1:length(v)
    end
    gs.x = gs.x[v_gs,:]
    gs.xi = gs.xi[v_gs,:]
    gs.weight = gs.weight[v_gs]
    gs.detJ = gs.detJ[v_gs]
    gs.Ni = gs.Ni[v_gs,:]
    gs.dNidx = gs.dNidx[v_gs,:]
    nlen = length(unique(element))
    mesh_picked = Mesh(type, dim, node, nlen, element, gs)
    return mesh_picked
end

function CombineMesh(meshes::Vector{Mesh})
    """
        combine meshes to a big mesh
        Inputs = meshes::Vector{Mesh}
        outputs = mesh_combined::Mesh
           e.g. meshnew = Combine([mesh1, mesh2, mesh3])
    """
    type = meshes[1].type
    dimension = meshes[1].dimension
    gsorder = meshes[1].gs.order
    for i = 2:length(meshes)
        if type != meshes[i].type
            error("the types of meshes do not match for combination!")
        end
        if dimension != meshes[i].dimension
            error("the dimensions of meshes do not match for combination!")
        end
        if gsorder != meshes[i].gs.order
            error("the orders of Gaussian quadrature do not match for combination!")
        end
    end

    n_mesh = length(meshes)
    n_element = zeros(Int64, n_mesh)
    n_node = zeros(Int64, n_mesh)
    n_gs = zeros(Int64, n_mesh)
    len_node = 0
    len_element = 0
    len_gs = 0
    for i = 1:n_mesh
        n_element[i] = size(meshes[i].element, 1)
        n_node[i] = meshes[i].nlen
        n_gs[i] = size(meshes[i].gs.x, 1)
        len_element += n_element[i]
        len_node += n_node[i]
        len_gs += n_gs[i]
    end

    element = zeros(Int64, len_element, size(meshes[1].element,2))
    node = zeros(Float64, len_node, size(meshes[1].node,2))
    x = zeros(Float64, len_gs, size(meshes[1].gs.x, 2))
    xi = zeros(Float64, len_gs, size(meshes[1].gs.xi, 2))
    weight = zeros(Float64, len_gs)
    detJ = zeros(Float64, len_gs)
    ele = zeros(Int64, len_gs)
    Ni = zeros(Float64, len_gs, size(meshes[1].gs.Ni, 2))
    dNidx = zeros(Float64, len_gs, size(meshes[1].gs.dNidx, 2))

    v_ele = 0
    v_node = 0
    v_gs = 0
    for i = 1:n_mesh
        element[v_ele + 1:v_ele + n_element[i],:] = meshes[i].element .+ v_node
        ele[v_gs + 1:v_gs + n_gs[i]]  = meshes[i].gs.ele .+ v_ele 
        node[v_node + 1:v_node + n_node[i],:] = meshes[i].node
        x[v_gs + 1:v_gs + n_gs[i], :] = meshes[i].gs.x
        xi[v_gs + 1:v_gs + n_gs[i], :]  = meshes[i].gs.xi
        weight[v_gs + 1:v_gs + n_gs[i]]  = meshes[i].gs.weight
        detJ[v_gs + 1:v_gs + n_gs[i]]  = meshes[i].gs.detJ
        Ni[v_gs + 1:v_gs + n_gs[i], :]  = meshes[i].gs.Ni
        dNidx[v_gs + 1:v_gs + n_gs[i], :]  = meshes[i].gs.dNidx
        v_ele += n_element[i]
        v_node += n_node[i]
        v_gs += n_gs[i]
    end

    gs = GaussPoint(x, xi, weight, detJ, ele, Ni, dNidx, gsorder)
    mesh_combined = Mesh(type, dimension, node, len_node, element, gs)
    return mesh_combined
end

function MultipleMesh(mesh::Mesh, n::Int64)
    """
        A function to duplicate a mesh by n times and combine them to one big mesh
            Inputs = mesh::Mesh
                    n::Int64
            Outputs = meshnew::Mesh
    """
    meshes = [mesh]
    for i = 1:n-1
        push!(meshes, mesh)
    end
    meshnew = CombineMesh(meshes)
    return meshnew
end

function GetGS(element::Array{Int64}, node::Array{Float64}, order::Int64, type::String, v=collect(1:size(element,1)))
    if type == "COH2D4"
        total_num = size(element, 1) * order
        x = zeros(Float64, total_num, 2)
        weight = zeros(Float64, total_num)
        detJ = zeros(Float64, total_num)
        ele = zeros(Int64, total_num)
        xi = zeros(Float64, total_num, 1)

        w, q = NCweight(order)
        count0 = 0
        for e = 1:size(element, 1)
            sctr = element[e, 1:4]

            x1 = node[sctr[1], 1]; y1 = node[sctr[1], 2]
            x2 = node[sctr[2], 1]; y2 = node[sctr[2], 2]
            x3 = node[sctr[3], 1]; y3 = node[sctr[3], 2]
            x4 = node[sctr[4], 1]; y4 = node[sctr[4], 2]

            x_m1 = 0.5 * (x1 + x4)
            y_m1 = 0.5 * (y1 + y4)
            x_m2 = 0.5 * (x2 + x3)
            y_m2 = 0.5 * (y2 + y3)

            dx = x_m2 - x_m1
            dy = y_m2 - y_m1
            L = sqrt(dx * dx + dy * dy)
            if L < 1e-15
                continue
            end

            for i = 1:length(w)
                count0 += 1
                xi_pt = q[i]
                N1 = 0.5 * (1.0 - xi_pt)
                N2 = 0.5 * (1.0 + xi_pt)
                x[count0, 1] = N1 * x_m1 + N2 * x_m2
                x[count0, 2] = N1 * y_m1 + N2 * y_m2
                weight[count0] = w[i]
                detJ[count0] = L / 2.0
                ele[count0] = v[e]
                xi[count0, 1] = xi_pt
            end
        end

        Ni, dNi = ShapeFunction2D(element, type, node, xi, ele)
        gs = GaussPoint(x, xi, weight, detJ, ele, Ni, dNi, order)
        return gs
    end

    if type == "L2" 
        dimen=1
        points = 1:2
    elseif type == "L3"
        dimen=1
        points = [1, 3]
    elseif type == "Q4"
        dimen=2
        points = 1:4
    elseif type == "B8"
        dimen=3
        points = 1:8
    end
    total_num = size(element,1) * order ^ dimen
    x = zeros(Float64, total_num ,dimen)
    weight = zeros(Float64, total_num)
    detJ = zeros(Float64, total_num) 
    ele = zeros(Int64, total_num)
    xi = zeros(Float64, total_num, dimen)
    w, q = GSweight(order,dimen)
    count0 = 0
    for e = 1:size(element, 1)
        sctr = element[e, points]
        for i = 1:size(w, 1)
            pt = q[i, :]
            N, dNdxi = LagrangeBasis(type, dimen, pt)
            J0 = dNdxi * node[sctr, 1:dimen]
            count0 = count0 + 1
            x[count0, 1:dimen] = N * node[sctr, 1:dimen]
            weight[count0] = w[i]
            detJ[count0] = det(J0)
            ele[count0] = v[e]
            xi[count0, 1:dimen] = pt
        end
    end
    if type in ["L2","L3"]
        Ni, dNi = ShapeFunction1D(element, type, node, xi, ele)
    elseif type == "Q4"
        Ni, dNi = ShapeFunction2D(element, type, node, xi, ele)
    else
        error("Unsupported element type $type for shape functions")
    end
    gs = GaussPoint(x, xi, weight, detJ, ele, Ni, dNi, order)
    return gs
end

function LagrangeBasis(type::String, dimen::Int64, coord::Array{Float64})
    N = zeros(Float64,2^dimen, 1)
    dNdxi = zeros(Float64, 2^dimen, dimen)
    if type == "L2" ||  type == "L3" 
        # 1------2 L2 TWO NODE LINE ELEMENT
        xi = coord[1]
        N[1,1] = (1.0 - xi)/2.0
        N[2,1] = (1.0 + xi)/2.0 
        dNdxi[1,1] = -1.0 /2.0  
        dNdxi[2,1] = 1.0 /2.0 
    elseif type == "Q4"
        ## Q4 FOUR NODE QURARILATERIAL ELEMENT
        # 4---3
        # |   |
        # 1---2
        xi = coord[1] 
        eta = coord[2]
        N=1/4 * [ (1-xi)*(1-eta);
            (1+xi)*(1-eta);
            (1+xi)*(1+eta);
            (1-xi)*(1+eta)]
        dNdxi=1/4 * [-(1-eta)   -(1-xi); 1-eta  -(1+xi); 1+eta  1+xi;   -(1+eta)    1-xi]
    elseif type == "B8"
        ## B4 EIGHT NODE BRICK ELEMENT
        # 4---3   8---7
        # |   |   |   |
        # 1---2 , 5---6
        xi=coord[1]
        eta=coord[2]
        zeta=coord[3]
        I1=1/2 - coord/2 
        I2=1/2 + coord/2
        N=[   I1[1]*I1[2]*I1[3];
            I2[1]*I1[2]*I1[3];
            I2[1]*I2[2]*I1[3];
            I1[1]*I2[2]*I1[3];
            I1[1]*I1[2]*I2[3];
            I2[1]*I1[2]*I2[3];
            I2[1]*I2[2]*I2[3];
            I1[1]*I2[2]*I2[3]   ]
        dNdxi=[-1+eta+zeta-eta*zeta   -1+xi+zeta-xi*zeta  -1+xi+eta-xi*eta;
            1-eta-zeta+eta*zeta   -1-xi+zeta+xi*zeta  -1-xi+eta+xi*eta;
            1+eta-zeta-eta*zeta    1+xi-zeta-xi*zeta  -1-xi-eta-xi*eta;
            -1-eta+zeta+eta*zeta    1-xi-zeta+xi*zeta  -1+xi-eta+xi*eta;
            -1+eta-zeta+eta*zeta   -1+xi-zeta+xi*zeta   1-xi-eta+xi*eta;
            1-eta+zeta-eta*zeta   -1-xi-zeta-xi*zeta   1+xi-eta-xi*eta;
            1+eta+zeta+eta*zeta   1+xi+zeta+xi*zeta   1+xi+eta+xi*eta;
            -1-eta-zeta-eta*zeta    1-xi+zeta-xi*zeta   1-xi+eta-xi*eta  ]/8
    end
    
    N=N'
    dNdxi=dNdxi'
    return N, dNdxi
end


function GSweight(order::Int64, dimen::Int64)
    if (order>10 || order<0)
        disp("Order of quadrature too high for Gaussian Quadrature")
    end
    r1pt = zeros(Float64, order) 
    r1wt = zeros(Float64, order)
    W = zeros(Float64, order^dimen)
    Q = zeros(Float64, order^dimen,dimen)
    if order == 1
        r1pt[1] = 0.000000000000000
        r1wt[1] = 2.000000000000000
        
    elseif order == 2
        r1pt[1] = 0.577350269189626
        r1pt[2] =-0.577350269189626
        
        r1wt[1] = 1.000000000000000
        r1wt[2] = 1.000000000000000
        
    elseif order ==  3
        r1pt[1] = 0.774596669241483
        r1pt[2] =-0.774596669241483
        r1pt[3] = 0.000000000000000
        
        r1wt[1] = 0.555555555555556
        r1wt[2] = 0.555555555555556
        r1wt[3] = 0.888888888888889
        
    elseif order ==  4
        r1pt[1] = 0.861134311594053
        r1pt[2] =-0.861134311594053
        r1pt[3] = 0.339981043584856
        r1pt[4] =-0.339981043584856
        
        r1wt[1] = 0.347854845137454
        r1wt[2] = 0.347854845137454
        r1wt[3] = 0.652145154862546
        r1wt[4] = 0.652145154862546
            
    elseif order ==  5
        r1pt[1] = 0.906179845938664
        r1pt[2] =-0.906179845938664
        r1pt[3] = 0.538469310105683
        r1pt[4] =-0.538469310105683
        r1pt[5] = 0.000000000000000
        
        r1wt[1] = 0.236926885056189
        r1wt[2] = 0.236926885056189
        r1wt[3] = 0.478628670499366
        r1wt[4] = 0.478628670499366
        r1wt[5] = 0.568888888888889
            
    elseif order ==  6
        r1pt[1] = 0.932469514203152
        r1pt[2] =-0.932469514203152
        r1pt[3] = 0.661209386466265
        r1pt[4] =-0.661209386466265
        r1pt[5] = 0.238619186003152
        r1pt[6] =-0.238619186003152
        
        r1wt[1] = 0.171324492379170
        r1wt[2] = 0.171324492379170
        r1wt[3] = 0.360761573048139
        r1wt[4] = 0.360761573048139
        r1wt[5] = 0.467913934572691
        r1wt[6] = 0.467913934572691
            
    elseif order ==  7
        r1pt[1] =  0.949107912342759
        r1pt[2] = -0.949107912342759
        r1pt[3] =  0.741531185599394
        r1pt[4] = -0.741531185599394
        r1pt[5] =  0.405845151377397
        r1pt[6] = -0.405845151377397
        r1pt[7] =  0.000000000000000
        
        r1wt[1] = 0.129484966168870
        r1wt[2] = 0.129484966168870
        r1wt[3] = 0.279705391489277
        r1wt[4] = 0.279705391489277
        r1wt[5] = 0.381830050505119
        r1wt[6] = 0.381830050505119
        r1wt[7] = 0.417959183673469
            
    elseif order ==  8
        r1pt[1] =  0.960289856497536
        r1pt[2] = -0.960289856497536
        r1pt[3] =  0.796666477413627
        r1pt[4] = -0.796666477413627
        r1pt[5] =  0.525532409916329
        r1pt[6] = -0.525532409916329
        r1pt[7] =  0.183434642495650
        r1pt[8] = -0.183434642495650
        
        r1wt[1] = 0.101228536290376
        r1wt[2] = 0.101228536290376
        r1wt[3] = 0.222381034453374
        r1wt[4] = 0.222381034453374
        r1wt[5] = 0.313706645877887
        r1wt[6] = 0.313706645877887
        r1wt[7] = 0.362683783378362
        r1wt[8] = 0.362683783378362
        
    else
        r1pt[1] = 0.9739065285
        r1pt[2] = -0.9739065285
        r1pt[3] =  0.8650633677
        r1pt[4] = -0.8650633677
        r1pt[5] =  0.6794095683
        r1pt[6] = -0.6794095683
        r1pt[7] = 0.4333953941
        r1pt[8] =-0.4333953941
        r1pt[9] = 0.1488743390
        r1pt[10] =-0.1488743390
        
        r1wt[1] =0.0666713443
        r1wt[2] =0.0666713443
        r1wt[3] =0.1494513492
        r1wt[4] =0.1494513492
        r1wt[5] =0.2190863625
        r1wt[6] =0.2190863625
        r1wt[7] = 0.2692667193
        r1wt[8] =0.2692667193
        r1wt[9] = 0.2955242247
        r1wt[10] = 0.2955242247
    end
    num=1
    if dimen == 1
        for i = 1:order
            Q[num, 1] = r1pt[i]
            W[num] = r1wt[i]
            num = num+1
        end
    elseif dimen == 2
        for i = 1:order
            for j = 1:order
                Q[num,1] = r1pt[i]   
                Q[num,2] = r1pt[j]
                W[num] = r1wt[i] * r1wt[j]
                num = num + 1
            end
        end
    else
        for i=1:order
            for j=1:order
                for k=1:order
                    Q[num,1] = r1pt[i]
                    Q[num,2] = r1pt[j]
                    Q[num,3] = r1pt[k]
                    W[num] = r1wt[i] * r1wt[j] * r1wt[k]
                    num = num + 1
                end
            end
        end
    end
    return W, Q
end

function NCweight(order::Int64)
    if order < 2 || order > 5
        error("Newton-Cotes order $(order) not supported (2-5)")
    end

    if order == 2
        q = [-1.0, 1.0]
        w = [1.0, 1.0]
    elseif order == 3
        q = [-1.0, 0.0, 1.0]
        w = [1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0]
    elseif order == 4
        q = [-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0]
        w = [1.0 / 4.0, 3.0 / 4.0, 3.0 / 4.0, 1.0 / 4.0]
    else
        q = [-1.0, -0.5, 0.0, 0.5, 1.0]
        w = [7.0 / 45.0, 32.0 / 45.0, 12.0 / 45.0, 32.0 / 45.0, 7.0 / 45.0]
    end
    return w, q
end

function ShapeFunction1D(element::Matrix{Int64}, type::String, node::Matrix{Float64}, xi::Array{Float64}, v::Vector{Int64})
    if type == "L3"
            # f1 = x-> (x .- 1).^2 / 4 
            # f2 =  x-> (1 .- x.^2) / 2
            # f3 =  x-> (x .+ 1).^2 / 4
            # df1 = x-> (x .- 1) / 2
            # df2 = x-> -x
            # df3 = x-> (x .+ 1) / 2
            ## another group of shape functions
            f1 = x-> 0.5 * x.^2 - 0.5 * x
            f2 =  x-> - x.^2 .+ 1
            f3 =  x-> 0.5 * x.^2 + 0.5 * x
            df1 = x-> x .- 1 / 2
            df2 = x-> -2 * x
            df3 = x-> x .+ 1 / 2
            
            ele_length = abs.(node[element[v, 3]] - node[element[v, 1]])
            Ni = cat(f1(xi), f2(xi), f3(xi),dims=2)
            dNidX = cat(df1(xi), df2(xi), df3(xi), dims=2)
            dXdx = 2 ./ ele_length * ones(1, 3)
            dNidx = dNidX .* dXdx
    elseif type == "L2"
            f1 =  x-> (1 .- x)/2 
            f2 =  x-> (1 .+ x)/2
            df1 =  x-> -0.5 * ones(Float64,size(x))
            df2 =  x-> 0.5 * ones(Float64,size(x))
            
            ele_length = abs.(node[element[v, 2]] - node[element[v, 1]])
            Ni = cat(f1(xi), f2(xi), dims=2)
            dNidX = cat(df1(xi), df2(xi), dims=2)
            dXdx = 2 ./ ele_length * ones(1, 2)
            dNidx = dNidX .* dXdx
    else
            error("Error: element type $(type) has not been implemented in ShapeFunction1D!\n")
    end
    return Ni, dNidx 
end

"""
    ShapeFunction2D: Q4 and COH2D4 shape functions and gradients.
    输入:
        element, type, node, xi (总高斯点×2), ele_map (高斯点对应的单元索引)
    输出:
        Ni: (ngs × 4)
        dNidx: (ngs × 8)  前4列 dN/dx, 后4列 dN/dy
"""
function ShapeFunction2D(element::Matrix{Int64}, type::String, node::Matrix{Float64}, xi::Array{Float64}, ele_map::Vector{Int64})
    if type == "Q4"
        total_gs = size(xi,1)
        nnode_ele = 4
        Ni = zeros(Float64, total_gs, nnode_ele)
        dNidx = zeros(Float64, total_gs, nnode_ele * 2)
        for g in 1:total_gs
            xi_g  = xi[g,1]; eta_g = xi[g,2]
            Nloc = 0.25 * [(1 - xi_g)*(1 - eta_g);
                           (1 + xi_g)*(1 - eta_g);
                           (1 + xi_g)*(1 + eta_g);
                           (1 - xi_g)*(1 + eta_g)]
            dN_dxi = 0.25 * [-(1 - eta_g)    -(1 - xi_g);
                              (1 - eta_g)    -(1 + xi_g);
                              (1 + eta_g)     (1 + xi_g);
                             -(1 + eta_g)     (1 - xi_g)]
            e = ele_map[g]
            sctr = element[e, :]
            x_e = node[sctr, 1:2]
            J = transpose(dN_dxi) * x_e
            invJ = inv(J)
            dN_dx = dN_dxi * invJ
            Ni[g, :] = vec(Nloc)'
            dNidx[g, 1:4] = dN_dx[:,1]'
            dNidx[g, 5:8] = dN_dx[:,2]'
        end
        return Ni, dNidx
    elseif type == "COH2D4"
        total_gs = size(xi,1)
        nnode_ele = 4
        Ni = zeros(Float64, total_gs, nnode_ele)
        dNidx = zeros(Float64, total_gs, nnode_ele * 2)
        for g in 1:total_gs
            xi_g = xi[g,1]
            N1 = 0.5 * (1.0 - xi_g)
            N2 = 0.5 * (1.0 + xi_g)

            e = ele_map[g]
            sctr = element[e, :]

            x1 = node[sctr[1], 1]; y1 = node[sctr[1], 2]
            x2 = node[sctr[2], 1]; y2 = node[sctr[2], 2]
            x3 = node[sctr[3], 1]; y3 = node[sctr[3], 2]
            x4 = node[sctr[4], 1]; y4 = node[sctr[4], 2]

            x_m1 = 0.5 * (x1 + x4)
            y_m1 = 0.5 * (y1 + y4)
            x_m2 = 0.5 * (x2 + x3)
            y_m2 = 0.5 * (y2 + y3)

            dx = x_m2 - x_m1
            dy = y_m2 - y_m1
            L = sqrt(dx * dx + dy * dy)
            if L < 1e-15
                continue
            end

            dN_dxi = [-0.5, 0.5]
            dN_ds = (2.0 / L) .* dN_dxi
            t_x = dx / L
            t_y = dy / L

            dN1_dx = dN_ds[1] * t_x
            dN1_dy = dN_ds[1] * t_y
            dN2_dx = dN_ds[2] * t_x
            dN2_dy = dN_ds[2] * t_y

            Ni[g, :] = [N1, N2, N2, N1]
            dNidx[g, 1:4] = [dN1_dx, dN2_dx, dN2_dx, dN1_dx]
            dNidx[g, 5:8] = [dN1_dy, dN2_dy, dN2_dy, dN1_dy]
        end
        return Ni, dNidx
    else
        error("Error: element type $type not implemented in ShapeFunction2D!")
    end
end
