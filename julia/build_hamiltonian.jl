include("./modular_arithmetic.jl")
include("./basis_construction.jl")

function get_momentum_projection(s::Int,nsites::Int,momentum::Int,base::Int)
    sᵢ = s
    smin = s
    R_to_min = 0

    R = 1 
    while R <= nsites
        sᵢ = circ_shift(sᵢ,1,;base=base,ndit=nsites)
        if sᵢ  < smin
            smin = sᵢ
            
            R_to_min = R 
        end 
        if sᵢ == s
            break 
        end
        R += 1
    end 
    # Check momentum compatibility
    (momentum*R) %nsites != 0 ? (-1,-1,-1) : (smin,R,R_to_min)
end

function apply_nn_hamiltonian(;H::Matrix{T}, site::Int,state_in::Int,nsites::Int, base::Int) where T <: Number
    hib_size = base^2 
    size(H) == (hib_size,hib_size) || throw("Hamiltonian shape $(size(H)) does not match that of the Hilbert space.")

    # Get the shifted_states
    state_shifted = circ_shift(state_in,-site;base=base,ndit=nsites)
    state_tail = mod(state_shifted, hib_size)
    state_head = state_shifted - state_tail 

    # Get the Hamiltonian coefficients 
    coeff = H[1:end,state_tail + 1] # Julia index start from 1
    
    # Get the shifted states 
    state_out = [circ_shift(state_head+i,site;base=base,ndit=nsites) for i in 0:(hib_size-1)]

    return coeff, state_out
end

# Hamiltonian without momentum 
struct NosymHamiltonianNN{T<: Number}
    Hnn::Matrix{T}
    nsites::Int
    base::Int
end

function (H::NosymHamiltonianNN{T})(v0::AbstractVector{T}) where {T <: Number}
    
    nsites = H.nsites
    base = H.base
    Hnn = H.Hnn
    Hv = zeros(T,base^nsites)

    for i in 1:length(v0)
        for site in 0:(nsites-1)
            state_in = i-1
            coeffs, bs = apply_nn_hamiltonian(;H=Hnn,site=site,state_in=state_in,nsites=nsites,base=base)
            for (c, b) in zip(coeffs, bs)
                Hv[b+1] += c*v0[i]
            end
        end
    end
    return Hv
end

block_size(H::NosymHamiltonianNN) = H.base^H.nsites
function dtype(H::NosymHamiltonianNN{T}) where {T}
    return T
end
function get_matrix(H::NosymHamiltonianNN{T}) where {T <: Number}
    dim = block_size(H)
    M = zeros(T,dim,dim)
    for i in 1:dim 
        e = zeros(T,dim)
        e[i] = one(T)
        M[1:end,i] = H(e)
    end
    return M
end

function get_nn_hamiltonian_nosym(hnn::Matrix{T},nsites::Int,base::Int) where {T <: Number}
    Id = Matrix{T}(I,base,base)
    H0 = hnn
    for _ in 1:(nsites-2)
        H0 = kron(H0,Id)
    end
    # Reshape into 2*nsites tensor 
    dims = ntuple(_-> base,2*nsites)
    H0 = reshape(H0,dims)

    # Build periodic sum 
    H_periodic = zero(H0)
    idcs = collect(1:nsites)

    for shift in 0:(nsites-1)

        left  = ((idcs .+ shift .- 1) .% nsites) .+ 1
        right = left .+ nsites

        perm = vcat(left, right)

        H_periodic .+= permutedims(H0, perm)
    end

    # ---- Reshape back to matrix ----
    dim = base^nsites
    return reshape(H_periodic, dim, dim)
end
# Periodic Hamiltonian 
struct PeriodicHamiltonianNN{T<: Number}
    Hnn::Matrix{T}
    momentum_sector::MomentumSector
end

function block_size(H::PeriodicHamiltonianNN)
    return length(H.momentum_sector.state_list)
end

function dtype(H::PeriodicHamiltonianNN{T}) where {T}
    return T
end

function (H::PeriodicHamiltonianNN{T})(v0::AbstractVector{T}) where {T <: Number}
    sector_size = block_size(H)
    length(v0) == sector_size || throw("v0 has shape $(size(v0)), but there are $(sector_size) basis")

    Hv = zeros(dtype(H),sector_size)
    basis = H.momentum_sector
    nsites = basis.nsites
    base = basis.base
    momentum = basis.momentum.k # Momentum (integer)
    charges = basis.charges
    state_list = basis.state_list
    period_list = basis.period_list

    k = 2*π*momentum/nsites
    for i in 1:length(v0) # Loop over columns 
        ψin = state_list[i]
        Pin = period_list[i]
        for site in 0:(nsites-1) # H(site, site+1) for site = 0, ..., site-1; site+nsites ∼ site 
            coeffs, ψouts = apply_nn_hamiltonian(;H=H.Hnn,site=site,state_in=ψin,nsites=nsites,base=base)
            # YM, 2/25/26. Possible optimization strategy: Shift just enoughto hit the representative, and check if representative exist, if exist, then just look up the period. You can probably also combine it with caching
            for (coeff,  ψout) in zip(coeffs, ψouts)
                ψout_rep, Pout,  ψout_r = get_momentum_projection(ψout,nsites,momentum,base)
                # println("R_to_min: $(ψout_r); k=$k")
                idc_out = find_index(ψout_rep,basis)
                if ψout_rep >= 0 && idc_out > 0
                    phase = exp(-im*k*ψout_r)
                    norm = sqrt(Pin/Pout)
                    # println("H[$idc_out] += $(norm*phase*coeff)")
                    Hv[idc_out] += norm*phase*coeff
                end
            end
        end
    end

    return Hv
end


let 
    using LinearAlgebra
    using KrylovKit
    Hnn = let 
        include("./spin_operators.jl")
        Sp,Sm,Sz = spin_operators(;S=3) # Spin 1 
        # (0.5+0*im)* (kron(Sp,Sm) + kron(Sm,Sp)) + kron(Sz,Sz)
        H0 = randn(ComplexF64,9,9)
        H0 + H0'
    end
    # Check that the basis constructed are the same as what we had in python 
    k0 = 2

    base = 3
    nsites = 12

    # get_nn_hamiltonian_nosym(Hnn,nsites,base)
    @time H = NosymHamiltonianNN(Hnn,nsites,base)
    N = block_size(H)
    println(N)
    # M = get_matrix(H)

    # @time vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),5,:SR)
    # vals2, vecs2 = eigen(M)

    # println("Using Krylov methods: ",real.(vals1))
    # println("Using matrix directly: ", vals2[1:10])

    # k_charge = Momentum(;k=k0,nsites=nsites,base=base)
    # z3_charge = ZNCharge(;N=3,c=0,nsites=nsites,charge_table=[0,1,-1])

    # sector = MomentumSector(;base=3,nsites=nsites,k=k0,charges=[z3_charge]);

    # H = PeriodicHamiltonianNN(Hnn,sector)
    # # print(H, dtype(H))
    # v0 = zeros(dtype(H),block_size(H))
    # @time H(v0)
end 