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
    function Momentum(;k::Int,nsites::Int,base::Int)
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
struct ZNCharge <: AbelianCharge
    N::Int
    c::Int
    nsites::Int
    charge_table::Vector{Int} # The basis element i has charge charge_table[i]

    function ZNCharge(;N::Int,c::Int,nsites::Int,charge_table::Vector{Int})
        length(charge_table) == N || throw("Charge table length should be the same as N")
        new(N,mod(c,N),nsites,[mod(cc,N) for cc in charge_table])
    end
end
Base.show(io::IO, c::ZNCharge) = print(io,"Z$(c.N)Charge(c=$(c.c),nsites=$(c.nsites),charge_table=$(c.charge_table))")


function get_charge(state::Int,charge::ZNCharge)
    str_state = num_to_digits(state,charge.N;ndigit=charge.nsites)
    tot_charge = 0
    for c in str_state
         tot_charge += charge.charge_table[c-'0'+1]
    end
    tot_charge % charge.N
    # sum([charge.charge_table[c-'0'+1] for c in str_state]) % charge.N # YM, 2/24 : Old version. One-liner but slightly slower.
end


has_charge(state::Int,charge::ZNCharge) = get_charge(state,charge) == charge.c

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

let 
    # Check that the basis constructed are the same as what we had in python 
    k0 = 0
    nsites = 15
    base = 3
    k_charge = Momentum(;k=k0,nsites=nsites,base=base)
    z3_charge = ZNCharge(;N=3,c=0,nsites=nsites,charge_table=[1,0,-1])

    @time basis = get_sector_basis(base,nsites,[k_charge,z3_charge])

end 