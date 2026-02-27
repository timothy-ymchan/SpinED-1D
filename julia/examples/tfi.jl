result = let 
    using LinearAlgebra
    using KrylovKit
    using Plots
    include("../build_hamiltonian.jl")
    include("../basis_construction.jl")

    # TFI model 
    function TFI(;g::Real)
        include("../spin_operators.jl")
        # Sp,Sm,Sz = spin_operators(ComplexF64;S=2) # Spin 1/2
        X = ComplexF64[0. 1.; 1. 0.]
        Z = ComplexF64[1. 0.; 0. -1.]
        Id = I(2)

        -kron(Z,Z) - g*kron(X,Id)
    end
    Hnn = TFI(;g=1)

    # Check that the basis constructed are the same as what we had in python 
    k0 = 2

    base = 2
    nsites = 8

    # get_nn_hamiltonian_nosym(Hnn,nsites,base)
    @time H = NosymHamiltonianNN(Hnn,nsites,base)
    N = block_size(H)
    println(N)
    # M = get_matrix(H)

    #println("Hermitian check: ", maximum(abs.(M-M')))
    num_eig = 20
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
        #M = get_matrix(H)
        # println("Hermitian? ", maximum(abs.(M .- M')))
        println("Sector size: ", block_size(H))
        println("Evaluation time: ")
        v0 = zeros(dtype(H),block_size(H))
        @time H(v0)
        
        # @time vals2, vecs2 = eigen(M)
        vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),5,:SR;ishermitian=true,tol=1e-12,krylovdim=20)
        push!(eigvals,real.(vals1))
    end
    eigvals = sort(vcat(eigvals...))
    scatter!(eigvals[1:num_eig],marker=:diamond,ms=3,label="Krylov (momentum)",color=:yellow)
    
    # Perform ground state scan 
    gs_energies = Dict()
    for nsites in 4:2:16
        gs_energy = []
        for g in 0:0.04:2

            eigvals = []
            for k in [0]
                Hnn = TFI(;g=g)
                sector = MomentumSector(;base=base,nsites=nsites,k=k,charges=Vector{AbelianCharge}());
                H = PeriodicHamiltonianNN(Hnn,sector)
                N = block_size(H)
                println("Sector size: ", block_size(H))
                println("Evaluation time: ")
                v0 = zeros(dtype(H),block_size(H))
                @time H(v0)
                vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),5,:SR;ishermitian=true,tol=1e-12,krylovdim=20)

                push!(eigvals,real.(vals1))
            end
            

            eigvals = sort(vcat(eigvals...))
            push!(gs_energy,real(minimum(eigvals)))
        end
        gs_energies[nsites] = gs_energy
    end
    gs_energies
end 

let 
    function derivatives(x,y)
        dx = diff(x)
        dy = diff(y)
        xbar = 0.5*(x[2:end].+x[1:end-1])
        xbar, dy ./dx
    end
    p = plot(xlabel="g",ylabel="-d2ϵ/dg2")
    for (k,v) in sort(collect(result); by = first)
        g = 0:0.04:2
        x1,y1 = derivatives(g,v./k)
        x2,y2 = derivatives(x1,y1)
        plot!(p,x2,-y2,label="nsites=$k",marker=:circ)
    end
    p
end