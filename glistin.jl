# code to process the SnowEx GLISTIN data

using StatsBase
using Clustering

"""Read GLISTIN annotation file."""
function readannot(filename::String)
    data = Dict()
    open(filename) do f
        for line in eachline(f)
            if !startswith(line, ';')
                tokens = split(line, ['(', '='])
                data[strip(tokens[1])] = strip(tokens[end])
            end
        end
    end
    return data
end

"""Read binary GLISTIN SCH data file."""
function loadfile(filename:: String)
    annot = readannot(replace(filename, "prc.sch", "ann"))
    nc = parse(Int, annot["SCH Number of Cross Track Samples"])
    nr = parse(Int, annot["SCH Number of Along Track Lines"])
    open(filename) do f
        data = read(f, Float32, (nc, nr))
    end
    return data'
end

"""Classify height error image block."""
function classify_block(h::Array{Float32, 2})
    nodata = -1.0
    trees = zeros(size(h))
    hmask0 = h .> nodata
    trees[.!hmask0] = nodata
    if any(hmask0)
        hf = percentile(h[hmask0], 95)
        hmaskf = h .< hf
        hh = h[hmask0 .& hmaskf]
        if length(hh) > 2
            k = kmeans(hcat(hh)', 2)
            classes = assignments(k)
            _, treeclass = findmax(k.centers) # tree class label is the largest in terms of error value
            ct = classes .== treeclass
            classes[:] = 0
            classes[ct] = 1
            trees[.!hmaskf] = 1
            trees[hmask0 .& hmaskf] = classes
        end
    end
    return trees .== 1
end

"""Classify height error image using a block-based K-Means
clustering algorithm."""
function classify_trees(error::Array{Float32, 2}, nx=300, ny=300, xincr=100, yincr=100)
    nr, nc = size(error)
    trees = zeros(Bool, size(error))
    for i in 1:yincr:nr
        for j in 1:xincr:nc
            i1 = i
            i2 = min(i+ny, nr)
            j1 = j
            j2 = min(j+nx, nc)
            h = error[i1:i2, j1:j2]
            trees[i1:i2, j1:j2] = classify_block(h)
            println("$i of $nr, $j of $nc")
        end
    end
    return trees
end
