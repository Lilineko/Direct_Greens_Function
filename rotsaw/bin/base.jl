rot = (x, r) -> mod1(x + r, 4)

if info["lattice type"] == 0
    include("bethe2.jl")
else
    include("square.jl")
end

struct Node
    x::Int64
    y::Int64
end

Path = Vector{Node}

D = Dict{Char, Node}(
    'E' => Node(1,0),
    'N' => Node(0,1),
    'W' => Node(-1,0),
    'S' => Node(0,-1)
)

invMove = Dict{String, String}(
    "E" => "W",
    "N" => "S",
    "W" => "E",
    "S" => "N"
)

directions = ["E", "N", "W", "S"]

iD = Dict{Int64, String}(
    0 => "E",
    1 => "N",
    2 => "W",
    3 => "S"
)

Di = Dict{String, Int64}(
    "E" => 0,
    "N" => 1,
    "W" => 2,
    "S" => 3
)

Graph = Vector{Vector{String}}

function addNode(a::Node, b::Node)
    return Node(a.x + b.x, a.y + b.y)
end

function castStateToPath(state::String = "")
    result = Path()
    push!(result, Node(0,0))
    for move in state
        push!(result, addNode(result[end], D[move]))
    end
    return result
end

function calculateGraph(graph::Graph)
    result = ones(ComplexF64, info["omega points"])
    t = info["tunneling constant"]

    doubles = 0.
    doSingle = false
    for it in 1:length(graph)
        vertex = graph[it]
        if vertex[1] == ""
            doSingle = true
            result ./= Gamma(vertex[1], vertex[2])
            continue
        end
        if doSingle
            result .*= t ./ Gamma(vertex[1], vertex[2])
        else
            parentVertex = graph[it + 1]
            moves = split(parentVertex[2], "")
            reachableMoves = reduce(*, [m for m in moves if m != vertex[1][end:end]])
            doubles = Gamma(parentVertex[1], reachableMoves) .- doubles .* t^2
            result .*= doubles ./ Gamma(vertex[1], vertex[2])
            doubles = 1. ./ doubles
        end
    end
    return result
end

function getOmegaRange()
    return [
        info["omega min"] + (info["omega max"] - info["omega min"]) * (it - 1) / (info["omega points"] - 1)
        for it in 1:info["omega points"]
    ]
end

function rotGF(m)
    G = zeros(ComplexF64, info["omega points"])
    if info["lattice type"] == 0

        bethe = getBetheWalks(length(m))
        for bra in bethe
            braPhase = conj(getPhase(bra, m))
            for ket in bethe
                ketPhase = getPhase(ket, m)
                G .+= (1 / length(bethe)) .* calculateGraph(getGraph(bra, ket)) .* (braPhase * ketPhase)
            end
        end
        # G = reduce(+,
        #     [
        #         calculateGraph(getGraph(key[1], key[2])) * value(m)
        #         for (key, value) in pairs(database[length(m)])
        #     ]
        # );

    else

        saws = getSelfAvoidingWalks(length(m))
        for bra in saws
            braPhase = conj(getPhase(bra, m))
            for ket in saws
                ketPhase = getPhase(ket, m)
                G .+= (1 / length(saws)) .* calculateGraph(getGraph(bra, ket)) .* (braPhase * ketPhase)
            end
        end

    end
    return G
end

function readData(fileLocation)
    rawData = JSON.parsefile(fileLocation, use_mmap = false)
    keyCollection = collect(keys(rawData))

    data = Dict{Vector{Int64}, Vector{ComplexF64}}()

    for key in keyCollection
        G = Vector{ComplexF64}(undef, length(getOmegaRange()))
        for it in 1:length(G)
            point = rawData[key][it]["re"] + im * rawData[key][it]["im"]
            G[it] = point
        end
        push!(data, eval(Meta.parse(key)) => G)
    end

    return data
end

nothing
