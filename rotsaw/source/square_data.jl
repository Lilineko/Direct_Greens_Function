factors = Dict{Int64, Vector{Function}}(
    0 => [
        m -> 1
    ]
    ,
    1 => [
        m -> (1 / 4) * reduce(+, [
            exp(-m[1] * im * 2pi * x / 4) *
            exp(m[1] * im * 2pi * x / 4)
            for
                x in 0:3
        ])
        ,
        m -> (1 / 4) * reduce(+, [
            exp(-m[1] * im * 2pi * x / 4) *
            exp(m[1] * im * 2pi * y / 4)
            for
                x in 0:3,
                y in 0:3
            if
                x != y
        ])
    ]
    ,
    2 => [
        m -> 1 / 3
        ,
        m -> (1 / 12) * reduce(+, [
            exp(-m[1] * im * 2pi * x / 4) *
            exp(-m[2] * im * 2pi * x / 3) *
            exp(m[1] * im * 2pi * x / 4) *
            exp(m[2] * im * 2pi * y / 3)
            +
            exp(-m[1] * im * 2pi * x / 4) *
            exp(-m[2] * im * 2pi * y / 3) *
            exp(m[1] * im * 2pi * x / 4) *
            exp(m[2] * im * 2pi * x / 3)
            for
                x in 0:3,
                y in 0:2
            if
                (y != x) &&
                (y != rot(x, 2))
        ])
        ,
        m -> (1 / 12) * reduce(+, [
            exp(-m[1] * im * 2pi * x / 4) *
            exp(-m[2] * im * 2pi * x / 3) *
            exp(m[1] * im * 2pi * y / 4) *
            exp(m[2] * im * 2pi * z / 3)
            +
            exp(-m[1] * im * 2pi * y / 4) *
            exp(-m[2] * im * 2pi * z / 3) *
            exp(m[1] * im * 2pi * x / 4) *
            exp(m[2] * im * 2pi * x / 3)
            for
                x in 1:4,
                y in 1:4,
                z in 1:3
            if
                (y != x) &&
                (y != rot(x, 2))
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
        ["EE", "NE"] => factors[2][3],
        ["EE", "NN"] => factors[2][4],
        ["EN", "EN"] => factors[2][5],
        ["EN", "ES"] => factors[2][6],
        ["EN", "NE"] => factors[2][7]
    )
)

nothing
