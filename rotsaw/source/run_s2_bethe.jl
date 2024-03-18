using JSON, OrderedCollections

info = OrderedDict{String, Union{Int64,Float64}}(
    "tunneling constant" => 1.0,
    "coupling constant" => 0.4,
    "broadening" => 0.01,
    "omega min" => -3.5,
    "omega max" => 12.5,
    "omega points" => 8001,
    "lattice type" => 0
)

include("base.jl")

magnons = 1000
JRange = [J for J in 0.01:0.005:1.0]
mRange = [Int64[], [0], [1], [0,0], [1,0], [1,1]]#,[0],[1],[2],[3],[0,0],[0,1],[1,0],[1,1],[2,0],[2,1],[3,1]]
isNoint = true

peaks = OrderedDict{Vector{Int64}, OrderedDict{Float64, Vector{Union{Vector{Float64}, Vector{ComplexF64}}}}}()
for m in mRange
    println("\n", " Calculating ", m, " in ", string(mRange))

    push!(peaks, m => OrderedDict{Float64, Vector{Vector{Float64}}}())
    for J in JRange
        println("\n", " >> Calculating J = ", J)

        dir = "./database_bethe/" * string(magnons) * "_" * string(length(m))
        if isNoint
            dir *= "_noint"
        end
        dir *= "/J=" * string(J) * "/"
        filename = string(m)
        G = readData(dir * filename)[m]

        # calculate G
        # info["coupling constant"] = J
        # G = rotGF(m)

        p, w = decompose(G)
        push!(peaks[m], J => [p, w, G])
    end
end


include("plotting.jl")
dir = "./figures/bethe/"
if isNoint
    dir *= "spc_noint/"
else
    dir *= "spc/"
end
plts = plotS(peaks)
for it in 1:length(plts)
    png(plts[it], dir * string(mRange[it]) * ".png")
    display(plts[it])
end
