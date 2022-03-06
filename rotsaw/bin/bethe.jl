function getOmegaRange()
    return [
        info["omega min"] + (info["omega max"] - info["omega min"]) * (it - 1) / (info["omega points"] - 1)
        for it in 1:info["omega points"]
    ]
end

function invG(state::String)
    J = info["coupling constant"]
    iδ = im * info["broadening"]

    path = castStateToPath(state)
    l = length(path)
    Ω = 2J * l #magnons and hole energy 2J for z = 4
    Ω -= J * max(0., l - 3/2) #interactions

    return getOmegaRange() .+ iδ .- Ω
end

function Sigma(state::String)
    t = info["tunneling constant"]
    J = info["coupling constant"]
    iδ = im * info["broadening"]
    maxDepth = 1000

    path = castStateToPath(state)
    l = length(path) + 1
    Ω = 2J * l
    Ω -= J * max(0., l - 3/2)
    Ω = getOmegaRange() .+ iδ .- Ω

    Ω .-= J * maxDepth
    result = Ω
    for _ in 1:maxDepth
        Ω .+= J
        result = Ω .- (3. * (t^2) ./ result)
    end
    result = (t^2) ./ result
end

function Gamma(k, s)
    invG(k) .- length(s) .* Sigma(k)
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
        if d == invMove[lastMove]
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

nothing
