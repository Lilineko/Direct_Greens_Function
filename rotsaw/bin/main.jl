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
include("plotting.jl")

JRange = [J for J in 0.0:0.005:1.0]
mRange = [[]]#[[0,0], [0,1], [1,0], [1,1]]
nMagnons = 10
overwrite = false

for J in JRange
    print("\n", ">> Calculating J = ", J)

    info["coupling constant"] = J
    dir = "./database/J=" * string(J)
    doJob = !isdir(dir)
    if doJob
        mkdir(dir)
    else
        doJob = overwrite
    end

    if doJob
        for m in mRange
            println("\n", ">>>> Calculating m = ", string(m))
            writeCpp(info)
            run(`g++ main2.cpp -I ../../_eigen_ -m64 -O3 -march=corei7-avx -std=c++17 -o main2.exe`)
            run(`./main2.exe $nMagnons`)

            data = OrderedDict{Vector{Int64}, Vector{ComplexF64}}()

            push!(data, m => rotGF(m))

            filename = string(m)
            file = open(dir * "/" * filename, "w")
            JSON.print(file, data, 2)
            close(file)
        end
    else
        println(" --- skipped")
    end
end

# data = OrderedDict{Vector{Int64}, Vector{ComplexF64}}()

# for x in 0:3
#     m = [x]
#     push!(data, m => rotGF(m))
# end

# for x in 0:1
#     for y in 0:1
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
