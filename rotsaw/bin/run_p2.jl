using JSON, OrderedCollections

info = OrderedDict{String, Union{Int64,Float64}}(
    "tunneling constant" => 1.0,
    "coupling constant" => 0.4,
    "broadening" => 0.01,
    "omega min" => -3.5,
    "omega max" => 12.5,
    "omega points" => 10001,
    "lattice type" => 1
)

include("base.jl")

nMagnons = 10
JRange = [J for J in 0.0:0.01:1.0]
mRange = [[0],[1]]#[[0,0], [0,1], [1,0], [1,1]]

peaks = OrderedDict{Vector{Int64}, OrderedDict{Float64, Vector{Vector{Float64}}}}()
for m in mRange
    println("\n", " Calculating ", m, " in ", string(mRange))

    push!(peaks, m => OrderedDict{Float64, Vector{Vector{Float64}}}())
    for J in JRange
        println("\n", " >> Calculating J = ", J)

        dir = "./database/" * string(nMagnons) * "/J=" * string(J) * "/"
        filename = string(m)
        G = readData(dir * filename)[m]

        p, w = decompose(G)
        push!(peaks[m], J => [p, w])
    end
end

# for x in 0:3
#     m = [x]
#     push!(data, m => rotGF(m))
# end

# for x in 0:3
#     for y in 0:2
#         m = [x, y]
#         push!(data, m => rotGF(m))
#     end
# end

# for x in 0:3
#     for y in 0:2
#         for z in 0:2
#             m = [x, y, z]
#             push!(data, m => rotGF(m))
#         end
#     end
# end

# file = open("./data/square.json", "w")
# JSON.print(file, data, 2)
# close(file)

# for (m, G) in pairs(data)
#     name = "[]"
#     if length(m) > 0
#         name = string(m)
#     end
#
#     dir = "./figures/square/"
#     if info["lattice type"] == 0
#         dir = "./figures/bethe/"
#     end
#
#     dir *= string(length(m)) * "/"
#     if !isdir(dir)
#         mkdir(dir)
#     end
#
#     png(plotA(G), dir * name * ".png")
# end

include("plotting.jl")
plts = plotPS2(peaks)
for it in 1:length(plts)
    png(plts[it], "./figures/" * string(mRange[it]) * ".png")
    display(plts[it])
end
