include("./modular_arithmetic.jl")
using DataFrames

struct PeriodTable
    table::DataFrame
    nsites::Int 
    base::Int
    function PeriodTable(table::Vector,nsites::Int,base::Int)
        new(DataFrame(table),nsites,base)
    end
    function PeriodTable(table::DataFrame,nsites::Int,base::Int)
        new(table,nsites,base)
    end
end

function get_period(s::Int, nsites::Int, base::Int = 3)
    R = 1
    si = s
    @debug "s0 $(to_base(si, base, nsites))"
    while R < nsites
        sip1 = circ_shift(si, 1, base, nsites)
        @debug "s$R $(to_base(sip1, base, nsites))"
        if sip1 < s
            return -1
        elseif sip1 == s
            return R
        end
        si = sip1
        R += 1
    end
    return R
end

function construct_period_table(nsites, base::Int =3)
    basis = []
    for s in 0:(base^nsites-1)
        P = get_period(s,nsites,base)
        if P > 0
            push!(basis,Dict(:basis=>s,:period=>P))
        end
    end
    sum(item[:period] for item in basis) == base^nsites || throw("The total number of basis does not add up to base^nsites = $(base^nsites)")
    return PeriodTable(basis,nsites,base)
end


function display(period_table::PeriodTable)
    print("Number of sites: $(period_table.nsites)\n")
    period_groups = groupby(period_table.table,:period)
    for subgrp in period_groups
        println("Period: $(subgrp[1,:period])")
        for rep in subgrp[:,:basis]
            println("\t$(to_base(rep,period_table.base,period_table.nsites))")
        end
    end
end

let
    period_table = construct_period_table(6)
    display(period_table)
end