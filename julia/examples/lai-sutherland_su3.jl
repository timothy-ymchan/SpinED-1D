let 
    using LinearAlgebra
    using KrylovKit
    using Plots
    using DataFrames
    using CSV

    include("../build_hamiltonian.jl")
    include("../basis_construction.jl")

    # Lai-Sutherland model with SU3 symmetry 
    Hnn = let 
        include("../spin_operators.jl")
        Sp,Sm,Sz = spin_operators(ComplexF64;S=3) # Spin 1
        heis = (0.5)* (kron(Sp,Sm) + kron(Sm,Sp)) + kron(Sz,Sz)
        # Lai-Sutherland model
        heis + heis*heis
    end

    base = 3
    nsites = 11

    # get_nn_hamiltonian_nosym(Hnn,nsites,base)
    @time H = NosymHamiltonianNN(Hnn,nsites,base)
    N = block_size(H)
    println(N)
    # M = get_matrix(H)

    #println("Hermitian check: ", maximum(abs.(M-M')))
    num_eig = 10
    # @time vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),num_eig,:SR;ishermitian=true,tol=1e-12,krylovdim=20)
    # @time vals2, vecs2 = eigen(M)

    # println("Using Krylov methods: ",real.(vals1))
    # println("Real and imaginary parts of Kyrlov: ", vals1)
    # println("Full matrix: ", real.(vals2)[1:num_eig])
    
    plot(xlabel="n",ylabel="Energy")
    # scatter!(real.(vals2)[1:num_eig],marker=:rect,ms=5,label="Full matrix")
    # scatter!(real.(vals1)[1:num_eig],marker=:circ,ms=5,label="Krylov (no sym)")

    # println("Using matrix directly: ", vals2[1:10])

    # Momentum and charge resolved
    eigvals = []
    U1a_charge_table = [1, -1, 0] # Charges for SU(3) basis 
    U1b_charge_table = [0, 1, -1] # Charges for SU(3) basis 
    
    function get_su3_weights(hw::Tuple{Int,Int})
        α1 = (2, -1)
        α2 = (-1, 2)
        weights = [hw]
        lw = [hw]

        while length(lw) > 0
            nw = []
            for w in lw 
                for k in 1:w[1]
                    push!(nw, (w[1] - k*α1[1], w[2] - k*α1[2]))
                end
                for k in 1:w[2]
                    push!(nw, (w[1] - k*α2[1], w[2] - k*α2[2]))
                end
            end
            nw = unique(nw)
            push!(weights, nw...)
            lw = nw
        end
        return unique(weights)
    end
    target_weights = weights = [
        get_su3_weights((0,0))...,  # 1
        get_su3_weights((1,0))...,  # 3
        get_su3_weights((0,1))...,  # 3'
        get_su3_weights((1,1))...,  # 8

        get_su3_weights((2,0))...,  # 6
        get_su3_weights((0,2))...,  # 6'
        get_su3_weights((2,1))...,  # 15a
        get_su3_weights((1,2))...,  # 15a
        get_su3_weights((2,2))...,  # 27

        
        get_su3_weights((3,0))...,  # 10
        get_su3_weights((0,3))...,  # 10'
        get_su3_weights((3,1))...,  # 24
        get_su3_weights((1,3))...,  # 24'
        # get_su3_weights((3,2))...,  # 42
        # get_su3_weights((2,3))...,  # 42'
        # get_su3_weights((3,3))...,  # 64

        get_su3_weights((4,0))...,  # 15b
        get_su3_weights((0,4))...,  # 15b
        # get_su3_weights((4,1))...,  # 35
        # get_su3_weights((1,4))...,  # 35'
        # get_su3_weights((4,2))...,  # 81
        # get_su3_weights((2,4))...,  # 81'
        # get_su3_weights((4,3))...,  # 154
        # get_su3_weights((3,4))...,  # 154'
        # get_su3_weights((4,4))...,  # 256

        get_su3_weights((5,0))...,  # 21
        get_su3_weights((0,5))...,  # 21
        # get_su3_weights((5,1))...,  # 48
        # get_su3_weights((1,5))...,  # 48'
        # get_su3_weights((5,2))...,  # 105
        # get_su3_weights((2,5))...,  # 105'
        # get_su3_weights((5,3))...,  # 210
        # get_su3_weights((3,5))...,  # 210'
        # get_su3_weights((5,4))...,  # 384
        # get_su3_weights((4,5))...,  # 384'
        # get_su3_weights((5,5))...,  # 441
    ]
    target_weights = unique(target_weights)

    # Filter out weights that is not possible given the number of sites 
    # The rule is p-q == L mod 3 for a chain of length L
    target_weights = filter(w -> (w[1] - w[2] - nsites) % 3 == 0, target_weights)
    Ntot = 0

    # Counting the number of states by purely enumerating all basis elements
    # U1a_charge = U1Charge(;c=0,nsites=nsites,charge_table=U1a_charge_table)
    # U1b_charge = U1Charge(;c=0,nsites=nsites,charge_table=U1b_charge_table)
    # weight_count = Dict()
    # for state in 0:(base^(nsites)-1)
    #     λ1 = get_charge(state,U1a_charge)
    #     λ2 = get_charge(state,U1b_charge)
    #     if (λ1,λ2) in keys(weight_count)
    #         weight_count[(λ1,λ2)] += 1
    #     else 
    #         weight_count[(λ1,λ2)] = 1
    #     end
    # end
    # println("Weight counts by enumeration: ")
    # weight_count = sort(collect(weight_count), by = x -> x[1])
    # for (weight, count) in weight_count
    #     println("\t$weight: $count")
    # end

    # Momentum and Cartan in SU(3)
    su3_weight_count = Dict(w => 0 for w in target_weights)
    #su3_weight_count = Dict(w => 0 for w in [-3,-2,-1,0,1,2,3])
    eigvals = []
    for k in 0:(nsites-1)
        for λ in target_weights

            U1a_charge = U1Charge(;c=λ[1],nsites=nsites,charge_table=U1a_charge_table)
            U1b_charge = U1Charge(;c=λ[2],nsites=nsites,charge_table=U1b_charge_table)
            sector = MomentumSector(;base=base, nsites=nsites,k=k,charges=[U1a_charge,U1b_charge]);
            H = PeriodicHamiltonianNN(Hnn,sector)
            N = block_size(H)
            if N > 0
                println("k = $k, λ = $λ sector [Size = $N] ")
            else
                continue 
            end
            Ntot += N
            su3_weight_count[λ] += N

            num_print = min(num_eig, N)
            if N < 2000
                # Diagonalize the matrix in full 
                M = get_matrix(H)
                println("\tHermitian? ", maximum(abs.(M .- M')))
                t = @elapsed vals2, vecs2 = eigen(M)
                println("\tFull diagonalization: ", real.(vals2)[1:num_print])
                println("\tTime elapsed: ", t)
                num_taken = min(num_eig,length(vals2))
                push!(eigvals,Dict("k"=>k,"U1a"=>λ[1],"U1b"=>λ[2],"vals"=>real.(vals2)[1:num_taken]))
            else
                # Use Krylov methods
                t= @elapsed vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),num_eig,:SR;ishermitian=true,tol=1e-12,krylovdim=50)
                println("\tKrylov diagonalization: ", real.(vals1)[1:num_print])
                println("\tTime elapsed: ", t)  
                num_taken = min(num_eig,length(vals1))
                push!(eigvals,Dict("k"=>k,"U1a"=>λ[1],"U1b"=>λ[2],"vals"=>real.(vals1)[1:num_taken]))
            end
        end
    end

    println("\nSummary:")
    println("Total size across all sectors: $Ntot")
    println("SU(3) weight counts: ")
    su3_weight_count = sort(collect(su3_weight_count), by = x -> x[1])
    for (weight, count) in su3_weight_count
        if count > 0 
            println("\t$weight: $count")
        end
    end

    df = DataFrame(eigvals)
    # df = DataFrame(eigvals)

    # # Saving the eigenvalues to a CSV file
    CSV.write("lai_sutherland_eigvals_su3-nsites-$(nsites).csv", df)
    
end 