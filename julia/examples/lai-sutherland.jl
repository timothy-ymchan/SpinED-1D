let 
    using LinearAlgebra
    using KrylovKit
    using Plots
    using DataFrames
    using CSV

    include("../build_hamiltonian.jl")
    include("../basis_construction.jl")

    # Lai-Sutherland model
    Hnn = let 
        include("../spin_operators.jl")
        Sp,Sm,Sz = spin_operators(ComplexF64;S=3) # Spin 1
        heis = (0.5)* (kron(Sp,Sm) + kron(Sm,Sp)) + kron(Sz,Sz)
        # Lai-Sutherland model
        heis + heis*heis
    end

    base = 3
    nsites = 15

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
    Sz_charge_table = [1,0,-1] # Charges for spin-1 basis 
    Sz_target = [0,-1,1,-2,2,-3,3,-4,4,-5,5] # Target Sz sectors to look at 
    for k in 0:(nsites-1)
        for Sz in Sz_target
            Sz_charge = U1Charge(;c=Sz,nsites=nsites,charge_table=Sz_charge_table)
            sector = MomentumSector(;base=base, nsites=nsites,k=k,charges=[Sz_charge]);
            H = PeriodicHamiltonianNN(Hnn,sector)
            N = block_size(H)
            println("k = $k, Sz = $Sz sector [Size = $N] ")
            if N < 2000
                # Diagonalize the matrix in full 
                M = get_matrix(H)
                println("\tHermitian? ", maximum(abs.(M .- M')))
                t = @elapsed vals2, vecs2 = eigen(M)
                println("\tFull diagonalization: ", real.(vals2)[1:5])
                println("\tTime elapsed: ", t)
                num_taken = min(num_eig,length(vals2))
                push!(eigvals,Dict("k"=>k,"Sz"=>Sz,"vals"=>real.(vals2)[1:num_taken]))
            else
                # Use Krylov methods
                t= @elapsed vals1, vecs1 = eigsolve(v->H(v),randn(dtype(H),N),num_eig,:SR;ishermitian=true,tol=1e-12,krylovdim=50)
                println("\tKrylov diagonalization: ", real.(vals1)[1:5])
                println("\tTime elapsed: ", t)  
                num_taken = min(num_eig,length(vals1))
                push!(eigvals,Dict("k"=>k,"Sz"=>Sz,"vals"=>real.(vals1)[1:num_taken]))
            end
        end
    end

    df = DataFrame(eigvals)

    # Saving the eigenvalues to a CSV file
    CSV.write("lai_sutherland_eigvals_$(nsites).csv", df)

    # Visualize the eigenvalues by momentum and Sz sectors 
    # p = plot(xlabel="k", ylabel="energy")
    # num_multiplet = 5 # Number of multiplets to show 
    # for ddf in groupby(df,[:k, :Sz])
    #     k0 = first(ddf.k)
    #     k0 = k0 > nsites ÷ 2 ? k0 - nsites : k0
    #     Sz = first(ddf.Sz)
    #     Sz_offset = 0.3 * Sz/ maximum(abs.(Sz_target)) # Offset points by Sz value for better visualization

    #     en = sort(vcat(ddf.vals...))
    #     tol = 1e-5 # tolerance for multiplets
    #     bins = []
    #     freq = []
    #     for i in 1:length(en)
    #         if i == 1
    #             push!(bins, en[i])
    #             push!(freq, 1)
    #         elseif abs(en[i] - bins[end]) < tol
    #             freq[end] += 1
    #         else
    #             push!(bins, en[i])
    #             push!(freq, 1)
    #         end
    #     end
        
    #     color = Sz == 0 ? :yellowgreen : (Sz > 0 ? :coral : :steelblue2)
    #     # marker = [:circ, :rect, :diamond, :utriangle, :dtriangle, :][abs.(Sz) + 1]
    #     marker = [:rect,:circ,  :diamond, :utriangle, :dtriangle,:cross, :xcross]

    #     ntake = min(num_multiplet,length(bins))

    #     markers = [marker[freq[i]] for i in 1:ntake]

    #     scatter!(fill(k0+Sz_offset,ntake),bins[1:ntake],marker=markers,ms=1.5,color=color,label=nothing)

    #     # Annotate degeneracies
    #     # for (bin, f) in zip(bins[1:num_multiplet], freq[1:num_multiplet])
    #     #     if f > 1
    #     #         annotate!(p, k0+Sz_offset, bin, text("x$f", 8 ,:red, :center))
    #     #     end
    #     # end
    # end
    # display(p)
    
end 