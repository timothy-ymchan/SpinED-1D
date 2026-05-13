let 
    using CSV 
    using Plots 
    using DataFrames 

    
    df = CSV.read("./lai_sutherland_eigvals_su3-nsites-12.csv", DataFrame)
    df = select!(df, [:k, :U1a, :U1b, :vals])
    df.vals = eval.(Meta.parse.(df[:, :vals]))

    # Filter out df with abs(U1a) > 3 or abs(U1b) > 3
    df = filter(row -> abs(row.U1a) <= 4 && abs(row.U1b) <= 4, df)

    nsites = maximum(df.k) + 1
    # Visualize the eigenvalues by momentum sectors and SU(3) weights
    p_list = []
    k_label = []
    num_multiplet = 5 # Number of multiplets to show
    for ddf in groupby(df,[:k])
        k0 = first(ddf.k)
        k0 = k0 > nsites ÷ 2 ? k0 - nsites : k0

        # Get the set of all multiplets
        weights = unique(ddf[:, [:U1a, :U1b]]) # Get unique SU(3) weights in this momentum sector
        weights = sort!(weights, [:U1a, :U1b]) # Sort by U1a and then U1b
        weight_labels = ["($(w.U1a),$(w.U1b))" for w in eachrow(weights)]
        println("k = $k0, SU(3) weights: ", weight_labels)

        p = plot(xlabel="weight (p,q)", ylabel="en [k=$k0]", xrotation = 45, title="L=$nsites")
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
        count = 0
        for (bin, f) in zip(bins[1:num_take], freq[1:num_take])
            println("\t$bin: [x$f]")
            hline!(p, [bin], color=:gray, linestyle=:dash, label=nothing)
            bin_str = string(round(bin, sigdigits=4))
            # xpos = count % 2 == 0 ? -1.5 : length(weight_labels) + 1.5 # Alternate annotation position to avoid overlap
            annotate!(-1.5, bin, text("en=$bin_str [x$f]", :grey, 5))
            count += 1
        end

        

        push!(p_list, p) 
        k_label = push!(k_label, k0)
    end

    # Saving results in a list of pdf files
    for (p, k) in zip(p_list, k_label)
        k_str = k >= 0 ? "k=+$k" : "k=$k"
        savefig(p, "./pdf_out/lai_sutherland_su3_$k_str.pdf")
    end

end