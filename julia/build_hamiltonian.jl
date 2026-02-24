include("./modular_arithmetic.jl")


function apply_nn_hamiltonian(;H::Matrix{T}, site::Int,state_in::Int,nsites::Int, base::Int) where T <: Number
    hib_size = base^2 
    size(H) == (hib_size,hib_size) || throw("Hamiltonian shape $(size(H)) does not match that of the Hilbert space.")

    # Get the shifted_states
    state_shifted = circ_shift(state_in,-site;base=base,ndit=nsites)
    state_tail = mod(state_shifted, hib_size)
    state_head = state_shifted - state_tail 

    # Get the Hamiltonian coefficients 
    coeff = H[:,state_tail + 1] # Julia index start from 1
    
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
end

function (H::PeriodicHamiltonianNN{T})(v0::AbstractVector{T}) where {T}
    Hv = zero(v0)
    for i in 1:length(v0) # Loop over columns 
        state_in = H.state_list[i]
        P_state_in = H.period_list[i]
        for site in 0:(H.nsites-1) # H(site, site+1) for site = 0, ..., site-1; site+nsites ∼ site 
            coeffs, state_outs = apply_nn_hamiltonian(;H=H.nn,site=site,state_in=state_in,nsites=H.nsites,base=H.base)
            for (coeff, state_out) in zip(coeffs,state_outs)
                
            end
        end
    end
end