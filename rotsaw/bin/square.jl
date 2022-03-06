function invG(state::String)
    J = info["coupling constant"]
    iδ = im * info["broadening"]

    path = castStateToPath(state)

    Ω = 2J * length(path) #magnons and hole energy 2J for z = 4
    for node in path
        for (_, move) in pairs(D)
            proximity = addNode(node, move)
            neighbours = reduce(+, [proximity == element for element in path[1:(end-1)]])
            # Ω -= 0.5J * neighbours #inlcude interactions -0.5J per interaction
        end
    end

    return getOmegaRange() .+ iδ .- Ω
end

function Gamma(k, s)
    result = invG(k)
    for move in s
        result .-= Sigma(k * move)
    end
    return result
end

function readSigma(state)
    J = info["coupling constant"]
    fileName = "output/J=" * rpad(string(J), 6, "0") * "/SE_" * state * ".txt"
    file = open(fileName, "r")
    lines = readlines(file)
    close(file)

    result = Vector{ComplexF64}(undef, length(lines))
    for it in 1:length(lines)
        line = lines[it]
        reim = split(split(split(line,")")[1],"(")[2],",")
        result[it] = parse(Float64, reim[1]) + im * parse(Float64, reim[2])
    end

    return result
end

function Sigma(state::String)
    moves = split(state, "")
    shift = findfirst(x -> x == moves[1], directions) - 1
    for it in 1:length(moves)
        position = findfirst(x -> x == moves[it], directions)
        moves[it] = directions[mod1(position - shift, 4)]
    end
    position = findfirst(x -> x != "E", moves)
    if !(position === nothing)
        if moves[position] == "S"
            for it in position:length(moves)
                if (moves[it] == "N") || (moves[it] == "S")
                    moves[it] = invMove[moves[it]]
                end
            end
        end
    end

    multiplicity = 1.
    if length(moves) == 1
        multiplicity = 4.
    else
        isAlmostStraight = reduce(&, moves[1:(end-1)] .== "E") && (moves[end] != "E")
        if isAlmostStraight
            multiplicity = 2.
        end
    end

    return readSigma(reduce(*, moves)) ./ multiplicity
end

function getReachableMoves(state, excludeSelf = false)
    if state == ""
        return reduce(*, directions)
    end
    result = ""
    for d in directions
        lastMove = state[end:end]
        if excludeSelf
            if d == lastMove
                continue
            end
        end
        path = castStateToPath(state)
        move = addNode(D[d[1]], path[end])
        if reduce(|, [move == node for node in path]) #if move d ends in path
            continue
        end
        result *= d
    end
    return result
end

function getGraph(bra::String, ket::String)
    result = Vector{Vector{String}}(undef, 0)
    minSize = min(length(bra), length(ket))

    l = 0
    for it in 1:minSize
        if bra[it:it] != ket[it:it]
            break
        end
        l += 1
    end

    commonStates = bra[1:l]
    leftStates = bra[(l+1):length(bra)]
    rightStates = ket[(l+1):length(ket)]


    state = ""
    vertex = [state, getReachableMoves(state)]
    push!(result, vertex)
    for it in 1:length(commonStates)
        state = commonStates[1:it]
        vertex = [state, getReachableMoves(state)]
        push!(result, vertex)
    end
    result = reverse(result)

    for it in 1:length(leftStates)
        state = commonStates * leftStates[1:it]
        vertex = [state, getReachableMoves(state)]
        push!(result, vertex)
    end

    for it in 1:length(rightStates)
        state = commonStates * rightStates[1:it]
        vertex = [state, getReachableMoves(state)]
        push!(result, vertex)
    end

    return result
end

function getSelfAvoidingWalks(nMoves)
    result = Vector{String}(undef, 0)
    function spanWalk(walk, depth)
        newWalk = deepcopy(walk)
        if depth != 0
            moves = split(getReachableMoves(newWalk), "")
            if moves[1] != ""
                for move in moves
                    spanWalk(newWalk * move, depth - 1)
                end
            end
        else
            push!(result, newWalk)
        end
    end
    spanWalk("", nMoves)
    return result
end

function getPhase(state::String, m)
    result = 1.
    nMoves = length(state)

    lastMove = 0
    lastState = ""
    for it in 1:nMoves
        moves = split(getReachableMoves(lastState), "")

        order = mod.([0, 1, 2, 3] .+ lastMove, 4)
        enum = [elem for elem in order if in(iD[elem], moves)]

        k = findfirst(x -> iD[x] == state[it:it], enum) - 1
        n = length(enum)

        result *= exp(-2pi * im * k * mod(m[it], n) / n)

        lastMove = Di[state[it:it]]
        lastState = state[1:it]
    end
    return result
end


function pick(x, y, xLim)
    indices = [i for (i, v) in enumerate(x) if v >= xLim[1] && v <= xLim[2]]
    return (x[indices], y[indices])
end

function decompose1(G)
    ω = getOmegaRange()
    A = real(1 ./ G)
    result = Vector{Float64}(undef, 0)
    isNegative = true
    for it in 2:length(ω)
        if isNegative
            if A[it] > 0
                isNegative = false
                push!(result, ω[it-1])
            end
        else
            if A[it] < 0
                isNegative = true
            end
        end
    end
    return result
end

function decompose2(G)
    ω = getOmegaRange()
    A = -imag(G) / pi
    result = Vector{Float64}(undef, 0)
    isRising = true
    for it in 2:length(ω)
        if isRising
            if A[it] - A[it-1] < 0
                isRising = false
                push!(result, ω[it-1])
            end
        else
            if A[it] - A[it-1] > 0
                isRising = true
            end
        end
    end
    return result
end

function decompose3(G)
    cutoff = 2.69 * 10^-7 / info["broadening"]
    ω = getOmegaRange()
    A = -imag(G) / pi
    result = Vector{Float64}(undef, 0)
    isRising = true
    for it in 2:length(ω)
        if isRising
            if A[it] - A[it-1] < 0
                if A[it-1] > cutoff
                    isRising = false
                    push!(result, ω[it-1])
                end
            end
        else
            if A[it] - A[it-1] > 0
                isRising = true
            end
        end
    end
    return result
end

function decompose4(G)
    cutoff = 10^-5 / info["broadening"]
    ω = getOmegaRange()
    A = -imag(G) / pi
    result = Vector{Float64}(undef, 0)
    isRising = true
    for it in 2:length(ω)
        if isRising
            if A[it] - A[it-1] < 0
                stipness = (2 * A[it-1] - A[it] - A[it-2]) / (ω[it] - ω[it-2])
                if stipness > cutoff
                    isRising = false
                    push!(result, ω[it-1])
                end
            end
        else
            if A[it] - A[it-1] > 0
                isRising = true
            end
        end
    end
    return result
end

function decompose5(G)
    cutoff = 10^-6 / info["broadening"]
    ω = getOmegaRange()
    A = -imag(G) / pi
    result = Vector{Float64}(undef, 0)
    isRising = true
    s, ωt = 0, 0
    for it in 2:length(ω)
        if isRising
            s += A[it] + A[it - 1]
            if A[it] - A[it-1] < 0
                if A[it-1] > cutoff
                    isRising = false
                    ωt = ω[it-1]
                end
            end
        else
            if A[it] - A[it-1] > 0
                isRising = true
                if s * (info["omega max"] - info["omega min"]) / (info["omega points"] - 1) > cutoff
                    push!(result, ωt)
                end
                s = 0
            end
            s += A[it] + A[it - 1]
        end
    end
    return result
end

function decompose(G)
    ω = getOmegaRange()
    A = -imag(G) / pi
    δ = info["broadening"]

    tt = Vector{Int64}(undef, 0)
    isRising = true
    for it in 2:length(ω)
        if isRising
            if A[it] - A[it-1] < 0
                isRising = false
                push!(tt, it-1)
            end
        else
            if A[it] - A[it-1] > 0
                isRising = true
            end
        end
    end

    splitMS = 1.0
    highPass = 10^-3 / (δ * pi)
    result = Vector{Float64}(undef, 0)
    weights = Vector{Float64}(undef, 0)
    for it in tt
        push!(result, ω[it])
        push!(weights, splitMS)
        if A[it] < highPass
            weights[end] *= (exp(1.0) * exp(-1.0 / (A[it] / highPass))) ^ 0.02
        end
    end

    return (result, weights)
end

function mapDecomposition(m, Jrange, info, nMagnons = 12)
    data = OrderedDict{Float64, Vector{Float64}}()
    for J in Jrange
        info["coupling constant"] = J
        println("\n", ">>> Calculating J = ", J)
        writeCpp(info)
        run(`g++ main2.cpp -I ../../_eigen_ -m64 -O3 -march=corei7-avx -std=c++17 -o main2.exe`)
        run(`./main2.exe $nMagnons`)

        G = rotGF(m)
        p = decompose(G)
        push!(data, J => p)
    end
    return data
end

function enum(peaks, m, J)
    data = peaks[m]
    p = data[J]
    ω = getOmegaRange()
    ip = [it for it in 1:length(ω) if ω[it] in p]
    return ip
end

nothing
