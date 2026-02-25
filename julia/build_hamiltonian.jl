include("./modular_arithmetic.jl")
include("./basis_construction.jl")

function get_momentum_projection(s::Int,nsites::Int,momentum::Int,base::Int)
    sᵢ = s
    smin = s
    R_to_min = 0
    for R in 1:nsites
        sᵢ = circ_shift(sᵢ,1,;base=base,ndit=nsites)
        if sᵢ  < smin
            smin = si 
            R_to_min = R 
        end 
        if sᵢ == s
            break 
        end
    end

    # Check momentum compatibility
    momentum %R %nsites != 0 ? (-1,-1,-1) : (smin,R,R_to_min)
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

struct PeriodicHamiltonianNN{T<: Number}
    Hnn::Matrix{T}
    base::Int
    nsites::Int
    momentum::Int
    state_list::Vector{Int}
    period_list::Vector{Int}
    charges::Vector{AbelianCharge}
    _index_cache::Dict{Int,Int}()
end

function PeriodicHamiltonianNN(Hnn::Matrix{T},momentum::Int,state_list::Vector{Int})
    
end

function find_index(H::PeriodicHamiltonianNN,state::Int)
    H._index_cache[state]
end

function (H::PeriodicHamiltonianNN{T})(v0::AbstractVector{T}) where {T}
    Hv = zero(v0)

    k = 2*π*H.momentum/H.nsites
    for i in 1:length(v0) # Loop over columns 
        ψin = H.state_list[i]
        Pin = H.period_list[i]
        for site in 0:(H.nsites-1) # H(site, site+1) for site = 0, ..., site-1; site+nsites ∼ site 
            coeffs, ψouts = apply_nn_hamiltonian(;H=H.nn,site=site,state_in=ψin,nsites=H.nsites,base=H.base)
            for (coeff,  ψout) in zip(coeffs, ψouts)
                ψout_rep, Pout,  ψout_r = get_momentum_projection(ψout,H.nsites,H.momentum,H.base)
                if ψout_rep >= 0 && all([has_charge(ψout,c) for c in H.charges])
                    ψout_rep in H.state_list || throw("The state $(ψout_rep) is missing from the list of states.")
                    phase = exp(-im*k*ψout_r)
                    norm = sqrt(Pin/Pout)
                    idc_out = find_index(H,ψout_rep)
                    Hv[idc_out] += norm*phase*coeff
                end
            end
        end
    end

    return Hv
end