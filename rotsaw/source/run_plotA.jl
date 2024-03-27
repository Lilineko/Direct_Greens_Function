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
include("plotting.jl")

nMagnons = 15
JRange = [J for J in 0.0:0.005:0.4]
mRange = [[0,0], [0,1], [1,0], [1,1]]

spc = OrderedDict{Vector{Int64}, OrderedDict{Float64, Vector{ComplexF64}}}()
for m in mRange
    println("\n", " Calculating ", m, " in ", string(mRange))

    push!(spc, m => OrderedDict{Float64, Vector{ComplexF64}}())
    for J in JRange
        println("\n", " >> Calculating J = ", J)

        dir = "./database/" * string(nMagnons) * "/J=" * string(J) * "/"
        filename = string(m)
        G = readData(dir * filename)[m]

        push!(spc[m], J => G)
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

for it in 1:length(mRange)
    for jt in 1:length(JRange)
        plt = plotA(spc[mRange[it]][JRange[jt]])
        png(plt, "./figures/square/spc/" * string(mRange[it]) * "_" * string(JRange[jt]) * ".png")
        display(plt)
    end
end
