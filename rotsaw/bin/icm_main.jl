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

struct Parameters
    nMagnons::Int64
    J::Float64
    mRange::Vector{Vector{Int64}}
end

function run_main(parameters::Parameters)
    nMagnons = parameters.nMagnons
    J = parameters.J
    mRange = parameters.mRange

    info["coupling constant"] = J
    dir = "./database/" * string(nMagnons) * "/J=" * string(J)
    if !isdir(dir)
        mkpath(dir)
    end

    outdir = "./output/J=" * rpad(string(J), 6, "0")
    if !isdir(outdir)
        mkpath(outdir)
    end

    run(`./main3.exe $nMagnons $J`)

    for m in mRange
        data = OrderedDict{Vector{Int64}, Vector{ComplexF64}}()
        push!(data, m => rotGF(m))

        filename = string(m)
        file = open(dir * "/" * filename, "w")
        JSON.print(file, data, 2)
        close(file)
    end
end

### System Parameters
nMagnons = 10
JRange = 0.0:0.1:1.0
mRange = [[]]

if length(ARGS) > 0
    nMagnons = eval(Meta.parse(ARGS[1]))
end
if length(ARGS) > 1
    JRange = eval(Meta.parse(ARGS[2]))
end
if length(ARGS) > 2
    mRange = eval(Meta.parse(ARGS[3]))
end

parameters = Vector{Parameters}()
for J in JRange
    push!(parameters, Parameters(nMagnons, J, mRange))
end

for it in 1:length(parameters)
    run_main(parameters[it])
end
