let 
    using CSV 
    using Plots 
    using DataFrames 

    
    df = CSV.read("./lai_sutherland_eigvals_su3-nsites-8.csv", DataFrame)
    df = select!(df, [:k, :U1a, :U1b, :vals])
    df.vals = eval.(Meta.parse.(df[:, :vals]))

    nsites = maximum(df.k) + 1
    # Visualize the eigenvalues by momentum sectors and SU(3) weights
    p_list = []
    num_multiplet = 5 # Number of multiplets to show
    for ddf in groupby(df,[:k])
        k0 = first(ddf.k)
        k0 = k0 > nsites ÷ 2 ? k0 - nsites : k0
   
        # Get the set of all multiplets
        weights = unique(ddf[:, [:U1a, :U1b]]) # Get unique SU(3) weights in this momentum sector
        weights = sort!(weights, [:U1a, :U1b]) # Sort by U1a and then U1b
        weight_labels = ["($(w.U1a),$(w.U1b))" for w in eachrow(weights)]
        println("k = $k0, SU(3) weights: ", weight_labels)

        p = plot(xlabel="weight (p,q)", ylabel="en [k=$k0]", xrotation = 45)
        p = xticks!(p, 1:length(weight_labels), weight_labels)

        en_all = []
        for wg in groupby(ddf, [:U1a, :U1b])

            U1a = first(wg.U1a)
            U1b = first(wg.U1b)
            en = sort(vcat(wg.vals...))

            # Binning the multiplets
            tol = 1e-5
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

            

            # Plotting 
            pos = findfirst(==("($U1a,$U1b)"), weight_labels)
            num_take = min(num_multiplet,length(bins))

            scatter!(fill(pos, num_take), bins[1:num_take], marker=:rect, ms=3,label=nothing, color=:steelblue2)

            # Annotate degeneracies
            for (bin, f) in zip(bins[1:num_take], freq[1:num_take])
                if f > 1
                    annotate!(pos, bin, text("x$f", :red, 8))
                end
            end

            push!(en_all, en)
        end
        en_all = sort(vcat(en_all...))
        println(en_all)

        # Summarize the multiplet structure in this momentum sector
        bins = []
        freq = []
        tol = 1e-5
        for i in 1:length(en_all)
            if i == 1
                push!(bins, en_all[i])
                push!(freq, 1)
            elseif abs(en_all[i] - bins[end]) < tol
                freq[end] += 1
            else
                push!(bins, en_all[i])
                push!(freq, 1)
            end
        end
        num_take = min(5,length(bins))
        println("Summary of multiplet structure for k=$k0:")
        for (bin, f) in zip(bins[1:num_take], freq[1:num_take])
            println("\t$bin: [x$f]")
            hline!(p, [bin], color=:gray, linestyle=:dash, label=nothing)
        end

        

        display(p)
        readline()
        # Get the set of all multiplets 
        
    end

end