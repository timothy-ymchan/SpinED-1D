let 
    using LinearAlgebra
    using KrylovKit
    using Plots
    include("../build_hamiltonian.jl")
    include("../basis_construction.jl")

    # Spin-1 Heisenberg model 
    Hnn = let 
        include("../spin_operators.jl")
        Sp,Sm,Sz = spin_operators(ComplexF64;S=3) # Spin 1/2
        (0.5)* (kron(Sp,Sm) + kron(Sm,Sp)) + kron(Sz,Sz)
    end

    # Check that the basis constructed are the same as what we had in python 
    k0 = 2

    base = 3
    nsites = 7

    # get_nn_hamiltonian_nosym(Hnn,nsites,base)
    @time H = NosymHamiltonianNN(Hnn,nsites,base)
    N = block_size(H)
    println(N)
    M = get_matrix(H)

    #println("Hermitian check: ", maximum(abs.(M-M')))
    num_eig = 10
    @time vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),num_eig,:SR;ishermitian=true,tol=1e-12,krylovdim=20)
    @time vals2, vecs2 = eigen(M)

    # println("Using Krylov methods: ",real.(vals1))
    # println("Real and imaginary parts of Kyrlov: ", vals1)
    # println("Full matrix: ", real.(vals2)[1:num_eig])
    
    plot(xlabel="n",ylabel="Energy")
    scatter!(real.(vals2)[1:num_eig],marker=:rect,ms=5,label="Full matrix")
    scatter!(real.(vals1)[1:num_eig],marker=:circ,ms=5,label="Krylov (no sym)")

    # println("Using matrix directly: ", vals2[1:10])

    eigvals = []
    for k in 0:(nsites-1)
        # k_charge = Momentum(;k=k,nsites=nsites,base=base)
        
        sector = MomentumSector(;base=base,nsites=nsites,k=k,charges=Vector{AbelianCharge}());
        H = PeriodicHamiltonianNN(Hnn,sector)
        N = block_size(H)
        #M = get_matrix(H)
        # println("Hermitian? ", maximum(abs.(M .- M')))
        println("Sector size: ", block_size(H))
        println("Evaluation time: ")
        v0 = zeros(dtype(H),block_size(H))
        @time H(v0)
        
        # @time vals2, vecs2 = eigen(M)
        vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),20,:SR;ishermitian=true,tol=1e-12,krylovdim=50)
        push!(eigvals,real.(vals1))
    end
    eigvals = sort(vcat(eigvals...))
    scatter!(eigvals[1:num_eig],marker=:diamond,ms=3,label="Krylov (momentum)",color=:yellow)
    
end 