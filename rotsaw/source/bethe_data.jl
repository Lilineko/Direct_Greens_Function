factors = Dict{Int64, Vector{Function}}(
    0 => [
        m -> 1
    ]
    ,
    1 => [
        m -> 1
        ,
        m -> 4 * (mod(m[1], 4) == 0) - 1
    ]
    ,
    2 => [
        m -> 1
        ,
        m -> (1 / 2) * (3 * (mod(m[2], 3) == 0) - 1)
        ,
        m -> (1 / 12) * reduce(+, [
            exp(im * (2pi / 4) * (y - v) * m[1]) *
            exp(im * (2pi / 3) * (z - x) * m[2])
            for
                v in 0:3,
                y in 0:3,
                x in 0:2,
                z in 0:2
            if
                (v != y) &&
                (x != rot(v, 2)) &&
                (z != rot(y, 2))
        ])
    ]
)

database = Dict{Int64, Dict}(
    0 => Dict{Vector{String}, Function}(
        ["", ""] => factors[0][1]
    )
    ,
    1 => Dict{Vector{String}, Function}(
        ["E", "E"] => factors[1][1],
        ["E", "N"] => factors[1][2]
    )
    ,
    2 => Dict{Vector{String}, Function}(
        ["EE", "EE"] => factors[2][1],
        ["EE", "EN"] => factors[2][2],
        ["EE", "NE"] => factors[2][3]
    )
)

nothing
