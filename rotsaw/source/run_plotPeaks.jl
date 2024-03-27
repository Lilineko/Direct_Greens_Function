using JSON, OrderedCollections

info = OrderedDict{String, Union{Int64,Float64}}(
    "tunneling constant" => 1.0,
    "coupling constant" => 0.4,
    "broadening" => 0.001,
    "omega min" => -3.5,
    "omega max" => 6.5,
    "omega points" => 10001,
    "lattice type" => 1
)

include("base.jl")
include("makeCpp.jl")

Jrange = [J for J in 0.4:0.01:0.4]

peaks = OrderedDict{Vector{Int64}, OrderedDict{Float64, Vector{Float64}}}()

mRange = [[0,0], [0,1], [1,0], [1,1]]
for m in mRange
    println("\n", "Calculating ", m, " in ", string(mRange))
    push!(peaks, m => mapDecomposition(m, Jrange, info, 10))
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
plt = plotP(peaks)
png(plt, "./figures/" * string(mRange) * ".png")
display(plt)
