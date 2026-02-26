using LinearAlgebra

function spin_operators(::Type{T};S::Int) where {T<: Number}
    # S = 2j+1; n = j-m where m = j,..,0,...,-j
    j = T(S-1)/2
    Sz = diagm([j-T(n) for n in 0:(S-1)])
    Sp,Sm = zero(Sz),zero(Sz)
    # <S,n1|S+|S,n2> = √(s-n2)n2 δ(n1,n2-1)
    # <S,n1|S-|S,n2> = √(s-n2-1)(n2+1) δ(n1,n2+1)
    for n1 in 0:S-2
        Sp[n1+1,n1+2] = sqrt((S-(n1+1))*(n1+1))
        Sm[n1+2,n1+1] = sqrt((S-(n1+1))*(n1+1))
    end
    Sp,Sm,Sz
end