# 1. Generate template for the structure of the Greens function to calculate
# 2. Run externall code to obtain needed data and upload it
# 3. Calculate the spectral function with rotations

using DelimitedFiles

info = Dict{String, Union{Int64,Float64}}(
    "tunneling constant" => 1.0,
    "coupling constant" => 0.4,
    "broadening" => 0.01,
    "omega min" => -6.0,
    "omega max" => 8.0,
    "omega points" => 1401
)

struct Gamma
    state::Vector{Char}
    moves::Vector{Char}
end

struct SEPattern
    nodes::Vector{Gamma}
    doubles::Vector{Gamma}
    singles::Vector{Float64}
end

struct GFCoefficient
    leftState::Vector{Char}
    rightState::Vector{Char}
    pattern::SEPattern
end

function getRepresentativePath(inPath::Vector{Char})::Vector{Char}
    D = ['E', 'N', 'W', 'S']
    path = deepcopy(inPath)
    if length(path) == 0
        return path
    end
    if path[1] != D[1]
        shift = findfirst(x -> x == path[1], D) - 1
        for jt in 1:length(path)
            position = findfirst(x -> x == path[jt], D)
            path[jt] = D[mod1(position - shift, 4)]
        end
    end
    it = 1
    while path[it] == D[1]
        it += 1
        if it > length(path)
            return path
        end
    end
    if path[it] == D[4]
        for jt = it:length(path)
            if path[jt] == D[2]
                path[jt] = D[4]
            elseif path[jt] == D[4]
                path[jt] = D[2]
            end
        end
    end
    return path
end

function invMove(move)
    D = ['E', 'N', 'W', 'S']
    return D[mod1(findfirst(x -> x == move, D) + 2, 4)]
end

function getStateList(nMoves::Int64)::Vector{Vector{Char}}
    D = ['E', 'N', 'W', 'S']
    stateList = [Char[]]
    if nMoves > 0
        stateList = Vector{Vector{Char}}(undef, 0)
        for move in D
            push!(stateList, [move])
        end
    end
    for it = 2:nMoves
        oldList = stateList
        stateList = Vector{Vector{Char}}(undef, 0)
        for state in oldList
            for move in D
                if invMove(move) != state[end]
                    temp = deepcopy(state)
                    push!(stateList, push!(temp, move))
                end
            end
        end
    end
    return stateList
end

function getSEPattern(leftState, rightState)::SEPattern
    D = ['E', 'N', 'W', 'S']
    nodes = Vector{Gamma}(undef, 0)
    doubles = Vector{Gamma}(undef, 0)
    singles = Vector{Float64}(undef, 0)
    maxLength = max(length(leftState), length(rightState))
    repLeftState = getRepresentativePath(leftState)
    repRightState = getRepresentativePath(rightState)
    push!(nodes, Gamma(Char[], D))
    for it in 1:maxLength
        lPass = it <= length(leftState)
        rPass = it <= length(rightState)
        if lPass && rPass
            isDouble = true
            state = repLeftState[1:it]
            moves = [move for move in D if move != invMove(state[end])]
            push!(nodes, Gamma(state, moves))
            if leftState[1:it] != rightState[1:it]
                isDouble = false
                state = repRightState[1:it]
                moves = [move for move in D if move != invMove(state[end])]
                push!(nodes, Gamma(state, moves))
            end
            if isDouble
                state = Char[]
                moves = [move for move in D if (move != repLeftState[it])]
                if length(doubles) > 0
                    state = repLeftState[1:it-1]
                    moves = [move for move in D if (move != repLeftState[it-1]) && (move != invMove(leftState[it]))]
                end
                push!(doubles, Gamma(state, moves))
            else
                push!(singles, 1)
                push!(singles, 1)
            end
        elseif lPass
            state = repLeftState[1:it]
            moves = [move for move in D if move != invMove(state[end])]
            push!(nodes, Gamma(state, moves))
            push!(singles, 1)
        elseif rPass
            state = repRightState[1:it]
            moves = [move for move in D if move != invMove(state[end])]
            push!(nodes, Gamma(state, moves))
            push!(singles, 1)
        end
    end
    return SEPattern(nodes, doubles, singles)
end

function getGFCoefficients(nMoves::Int64)::Vector{GFCoefficient}
    result = Vector{GFCoefficient}(undef, 0)
    stateList = getStateList(nMoves)
    for leftState in stateList
        for rightState in stateList
            pattern = getSEPattern(leftState, rightState)
            push!(result, GFCoefficient(leftState, rightState, pattern))
        end
    end
    return result
end

function getInputData(GFCoefficients::Vector{GFCoefficient})
    SEList = Set(Vector{Char}[])
    for c in GFCoefficients
        for node in c.pattern.nodes
            for move in node.moves
                path = deepcopy(node.state)
                if length(path) == 0
                    SEPath = append!(path, 'E')
                    push!(SEList, SEPath)
                else
                    isStraight = true
                    for element in path
                        isStraight &= element == 'E'
                    end
                    if isStraight
                        if move == 'S'
                            move = 'N'
                        end
                    end
                    SEPath = append!(path, move)
                    push!(SEList, SEPath)
                end
            end
        end
        for double in c.pattern.doubles
            for move in double.moves
                path = deepcopy(double.state)
                if length(path) == 0
                    SEPath = append!(path, 'E')
                    push!(SEList, SEPath)
                else
                    isStraight = true
                    for element in path
                        isStraight &= element == 'E'
                    end
                    if isStraight
                        if move == 'S'
                            move = 'N'
                        end
                    end
                    SEPath = append!(path, move)
                    push!(SEList, SEPath)
                end
            end
        end
    end
    return SEList
end

function makeInputFile(inputData, fileLocation = "./")
    file = open(fileLocation * "input.txt", "w")
    for id in inputData
        write(file, reduce(*, id)) + write(file, "\n")
    end
    close(file)
end

function runSECalculations(inputData, nMagnons)
    # check if files exists and only calculate if not
    makeInputFile(inputData)
    run(`./main.exe $nMagnons`)
end

function calculateSEBethe(inputData, fileLocation = "./output_bethe/")
    for id in inputData
        l = length(id)
        fileName = fileLocation * "SE_" * reduce(*, id) * ".txt"
        ω = getOmegaRange()
        iδ = im * info["broadening"]
        J = info["coupling constant"]
        t = info["tunneling constant"]
        Ω = ω .+ iδ .- J * (2 + l)
        if l > 0
            Ω .-= 0.5 * J
        end
        maxDepth = 1000
        Ω .-= J * maxDepth
        result = Ω
        for _ in 1:maxDepth
            Ω .+= J
            result = Ω .- (3. * t * t ./ result)
        end
        result = 1. ./ result
        writable = ["(" * string(real(v)) * "," * string(imag(v)) * ")" for v in result]
        file = open(fileName, "w")
        writedlm(file, writable, '\n')
        close(file)
    end
end

function readSEData(inputData, fileLocation = "./")::Dict{Vector{Char}, Vector{ComplexF64}}
    result = Dict{Vector{Char}, Vector{ComplexF64}}()
    for id in inputData
        fileName = fileLocation * "SE_" * reduce(*, id) * ".txt"
        file = open(fileName, "r")
        lines = readlines(file)
        close(file)
        se = Vector{ComplexF64}(undef, length(lines))
        for it in 1:length(lines)
            line = lines[it]
            reim = split(split(split(line,")")[1],"(")[2],",")
            se[it] = parse(Float64, reim[1]) + im * parse(Float64, reim[2])
        end
        push!(result, id => se)
    end
    return result
end

function getOmegaRange()
    return [
        info["omega min"] + (info["omega max"]-info["omega min"]) * (it-1) / (info["omega points"]-1)
        for it in 1:info["omega points"]
    ]
end

function getInvBareGF(Γ::Gamma)
    ω = getOmegaRange()
    iδ = im * info["broadening"]
    J = info["coupling constant"]
    Ω = ω .+ iδ .- J * (2 + length(Γ.state))
    if length(Γ.state) > 0
        Ω .-= 0.5 * J
    end
    return Ω
end

function buildGamma(Γ::Gamma, SEData)
    if length(Γ.state) == 0
        return (getInvBareGF(Γ) .- (length(Γ.moves) .* SEData[['E']]))
    end
    result = getInvBareGF(Γ)
    for move in Γ.moves
        isStraight = true
        for element in Γ.state
            isStraight &= element == 'E'
        end
        if isStraight
            if move == 'S'
                move = 'N'
            end
        end
        path = deepcopy(Γ.state)
        append!(path, move)
        result .-= SEData[path]
    end
    return result
end

function calculateGFCoefficient(GFCoefficient, SEData)
    c = GFCoefficient
    t = info["tunneling constant"]
    #loop over nodes
    nodesValue = 1.
    for node in c.pattern.nodes
        nodesValue = buildGamma(node, SEData) .* nodesValue
    end
    #loop over doubles
    doublesValue = 1.
    CF = 0.
    for double in reverse(c.pattern.doubles)
        # continued fraction
        Γ = buildGamma(double, SEData)
        CF = Γ .- (t^2) .* CF
        doublesValue = CF .* doublesValue
        CF = 1. ./ CF
    end
    singlesValue = 1.
    singlesValue = t ^ length(c.pattern.singles)
    return (doublesValue .* singlesValue) ./ nodesValue
end

function getGFList(gfc, SEData)
    D = ['E', 'N', 'W', 'S']
    ### dict with left, right state, gfc value, and phase factor/s
    result = Dict{Vector{Vector{Char}}, Vector{ComplexF64}}()
    for it in 1:length(gfc)
        key = [gfc[it].leftState, gfc[it].rightState]
        value = calculateGFCoefficient(gfc[it], SEData)
        push!(result, key => value)
    end
    return result
end

function getPhase(state, mList)
    D = ['E', 'N', 'W', 'S']
    phase = 1.
    if length(state) > 0
        path = deepcopy(state)
        dist = findfirst(x -> x == path[1], D) - 1
        phase *= exp(-im * 0.5 * pi * dist * mList[1])
        for it = 2:length(path)
            # lastMove = path[it - 1]
            # move = path[it]
            # offset = findfirst(x -> x == lastMove, D)
            # dist = 0
            # for jt in offset:(offset + 3)
            #     index = mod1(jt, 4)
            #     if D[index] != invMove(lastMove)
            #         if D[index] == move
            #             break
            #         end
            #         dist += 1
            #     end
            # end
            dist = findfirst(x -> x == path[it], D) - 1
            phase *= exp(-im * 2pi * (1. / 3.) * dist * mList[it])
        end
    end
    return phase
end

function calculateRotGF(mList)
    nMoves = length(mList)
    gfc = getGFCoefficients(nMoves)
    inputData = getInputData(gfc)
    calculateSEBethe(inputData)
    SEData = readSEData(inputData, "output_bethe/")
    GFList = getGFList(gfc, SEData)
    result = zeros(ComplexF64, length(getOmegaRange()))
    for (key, value) in pairs(GFList)
        ls, rs = key[1], key[2]
        lph, rph = getPhase(ls, mList), getPhase(rs, mList)
        result .+= conj(lph) * rph .* value
    end
    return result
end
