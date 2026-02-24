include("./modular_arithmetic.jl")

function get_period(s::Int,nsites::Int;base::Int)
    R = 1 # Period 
    sᵢ = s
    while R < nsites
        sᵢ₊₁ = circ_shift(sᵢ,1;base=base,ndit=nsites)
        if sᵢ₊₁ < s 
            return -1
        elseif sᵢ₊₁ == s 
            return R 
        end
        sᵢ = sᵢ₊₁
        R += 1
    end
    return R
end

abstract type AbelianCharge end

# Crystal momentum 
struct Momentum <: AbelianCharge 
    nsites::Int
    k::Int
    base::Int
    function Momentum(k::Int,nsites::Int;base::Int)
        new(nsites,mod(k,nsites),base)
    end
end

Base.show(io::IO, m::Momentum)=print(io,"Momentum(k=$(m.k),nsites=$(m.nsites);base=$(m.base))")

function has_charge(state::Int,charge::Momentum)
    P = get_period(state,charge.nsites;base=charge.base)
    if P < 0
        return false
    end 
    return charge.k * P % charge.nsites == 0
end

# Zn charge 
# struct ZnCharge <: AbelianCharge
#     n::Int
#     charge_table::Vector{}
# end

# Compute sector basis 
function get_sector_basis(base::Int,nsites::Int,charges::Vector{<:AbelianCharge})
    state_list = Vector{Int}()
    for state in 0:(base^(nsites)-1)
        if all([has_charge(state,c) for c in charges])
            push!(state_list,state)
        end
    end
    return sort(state_list)
end