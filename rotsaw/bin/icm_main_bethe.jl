using JSON, OrderedCollections

info = OrderedDict{String, Union{Int64,Float64}}(
    "tunneling constant" => 1.0,
    "coupling constant" => 0.4,
    "broadening" => 0.02,
    "omega min" => -3.5,
    "omega max" => 12.5,
    "omega points" => 4001,
    "lattice type" => 0
)

include("base.jl")

struct Parameters
    nMagnons::Int64
    J::Float64
    mRange::Vector{Vector{Int64}}
end

function run_main_bethe(parameters::Parameters)
    nMagnons = parameters.nMagnons
    J = parameters.J
    mRange = parameters.mRange

    global maxDepth = min(1000, round(10 / J))

    info["coupling constant"] = J
    println("\n", " >> Calculating J = ", J)

    for m in mRange
        data = OrderedDict{Vector{Int64}, Vector{ComplexF64}}()
        push!(data, m => rotGF(m))

        dir = "./database_bethe/" * string(nMagnons) * "_" * string(length(m))
        if isNoint
            dir *= "_noint"
        end
        dir *= "/J=" * string(J)
        if !isdir(dir)
            mkpath(dir)
        end

        filename = string(m)
        file = open(dir * "/" * filename, "w")
        JSON.print(file, data, 2)
        close(file)
    end
end

### System Parameters
maxDepth = 1000 ### cmp `global maxDepth` in `run_main_bethe`
JRange = 0.01:0.005:1.0
mRange = [Int64[],[0],[1],[0,0],[0,1],[1,1]]
isNoint = false

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
    push!(parameters, Parameters(maxDepth, J, mRange))
end

for it in 1:length(parameters)
    run_main_bethe(parameters[it])
end
