let 
    using LinearAlgebra
    using KrylovKit
    using Plots
    include("../build_hamiltonian.jl")
    include("../basis_construction.jl")

    # Spin-1 Heisenberg model 
    Hnn = let 
        include("../spin_operators.jl")
        Sp,Sm,Sz = spin_operators(ComplexF64;S=3) # Spin 1
        heis = (0.5)* (kron(Sp,Sm) + kron(Sm,Sp)) + kron(Sz,Sz)
        # Lai-Sutherland model
        heis + heis*heis
    end

    base = 3
    nsites = 9

    # get_nn_hamiltonian_nosym(Hnn,nsites,base)
    @time H = NosymHamiltonianNN(Hnn,nsites,base)
    N = block_size(H)
    println(N)
    # M = get_matrix(H)

    #println("Hermitian check: ", maximum(abs.(M-M')))
    num_eig = 30
    # @time vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),num_eig,:SR;ishermitian=true,tol=1e-12,krylovdim=20)
    # @time vals2, vecs2 = eigen(M)

    # println("Using Krylov methods: ",real.(vals1))
    # println("Real and imaginary parts of Kyrlov: ", vals1)
    # println("Full matrix: ", real.(vals2)[1:num_eig])
    
    plot(xlabel="n",ylabel="Energy")
    # scatter!(real.(vals2)[1:num_eig],marker=:rect,ms=5,label="Full matrix")
    # scatter!(real.(vals1)[1:num_eig],marker=:circ,ms=5,label="Krylov (no sym)")

    # println("Using matrix directly: ", vals2[1:10])

    eigvals = []
    for k in 0:(nsites-1)
        # k_charge = Momentum(;k=k,nsites=nsites,base=base)
        
        sector = MomentumSector(;base=base,nsites=nsites,k=k,charges=Vector{AbelianCharge}());
        H = PeriodicHamiltonianNN(Hnn,sector)
        N = block_size(H)
        # M = get_matrix(H)
        # println("Hermitian? ", maximum(abs.(M .- M')))
        println("Sector size: ", block_size(H))
        println("Evaluation time: ")
        v0 = zeros(dtype(H),block_size(H))
        @time H(v0)
        
        # @time vals2, vecs2 = eigen(M)
        println("Krylov time: ")
        @time vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),20,:SR;ishermitian=true,tol=1e-12,krylovdim=50)
        push!(eigvals,real.(vals1))
    end

    # Visualize the eigenvalues by momentum sectors
    p = plot(xlabel="k", ylabel="energy")
    for k in 0:(nsites-1)
        k0 = k > nsites ÷ 2 ? k - nsites : k
        en = sort(eigvals[k+1])
        

        tol = 1e-5 # tolerance for multiplets 
        bins = []
        freq = []
        
        for i in 1:length(en)
            if i == 1
                push!(bins, en[i])
                push!(freq, 1)
            elseif abs(en[i] - bins[end]) < tol
                freq[end] += 1
            else
                push!(bins, en[i])
                push!(freq, 1)
            end
        end

        kk0 = fill(k0, length(bins))
        scatter!(p, kk0, bins, marker=:rect, ms=2, label=nothing, color=:royalblue)
        # Plot degeneracies as annotations 
        for (bin, f) in zip(bins, freq)
            if f > 1
                annotate!(p, k0, bin, text("x$f", 8 ,:red, :center))
            end
        end
        # scatter!(p,kk0, en, marker=:rect, ms=2,label=nothing,color=:red)
        
    end
    display(p)
    #eigvals = sort(vcat(eigvals...))
    #scatter!(eigvals[1:num_eig],marker=:diamond,ms=3,label="Krylov (momentum)",color=:yellow)
    
end 