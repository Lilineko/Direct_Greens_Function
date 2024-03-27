## approach based on definitions
using DelimitedFiles

info = Dict{String, Union{Int64,Float64}}(
    "tunneling constant" => 1.0,
    "coupling constant" => 0.4,
    "broadening" => 0.05,
    "omega min" => -6.0,
    "omega max" => 8.0,
    "omega points" => 1401
)

function getOmegaRange()
    return [
        info["omega min"] + (info["omega max"]-info["omega min"]) * (it-1) / (info["omega points"]-1)
        for it in 1:info["omega points"]
    ]
end

function Gamma(k, s)
    invG(k) - s * Sigma(k)
end

function invG(k)
    ω = getOmegaRange()
    iδ = im * info["broadening"]
    J = info["coupling constant"]
    Ω = ω .+ iδ .- J * (2 + k)
    if k > 0
        Ω .-= 0.5 * J
    end
    return Ω
end

function Sigma(k)
    iδ = im * info["broadening"]
    J = info["coupling constant"]
    t = info["tunneling constant"]
    Ω = invG(k + 1)
    maxDepth = 1000
    Ω .-= J * maxDepth
    result = Ω
    for _ in 1:maxDepth
        Ω .+= J
        result = Ω .- (3. * (t^2) ./ result)
    end
    result = (t^2) ./ result
end

G = 1. ./ Gamma(0, 4)

t = info["tunneling constant"]

G_EE_EE = (Gamma(0, 3) .- ((t^2) ./ Gamma(1, 2))) .* Gamma(1, 2) ./ (Gamma(0, 4) .* Gamma(1, 3) .* Gamma(2, 3))

# G_EE_EE = (Gamma(0, 3) .* Gamma(1, 2) .- t^2) ./ (Gamma(0, 4) .* Gamma(1, 3) .* Gamma(2, 3))

G_EE_EN = (Gamma(0, 3) .* t^2) ./ (Gamma(0, 4) .* Gamma(1, 3) .* Gamma(2, 3) .* Gamma(2, 3))

G_EE_NE = (t^4) ./ (Gamma(0, 4) .* Gamma(1, 3) .* Gamma(1, 3) .* Gamma(2, 3) .* Gamma(2, 3))

mList = [0,0]
fc3 = (1/12)*sum([exp(im * 0.5pi * (y - v) * mList[1]) * exp(im * 2pi * (1/3) * (z - x) * mList[2]) for v in 1:4, y in 1:4, x in 1:3, z in 1:3 if ((v != y) && (x != mod1(v + 2, 4)) && (z != mod1(y + 2, 4)))])
G = (1.) .* G_EE_EE .+ (3/2 * (mList[2]==0) - 1/2) .* G_EE_EN .+ fc3 .* G_EE_NE

# display(plot(getOmegaRange(), -imag(G)/pi))

display(plot(getOmegaRange(), -imag(G)/pi))
