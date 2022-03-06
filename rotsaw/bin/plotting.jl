using Plots, Glob, Plots.PlotMeasures, LaTeXStrings

colors = Dict{Vector{Int64}, String}(
    [0] => "blue",
    [1] => "orange",
    [0,0] => "blue",
    [0,1] => "green",
    [1,0] => "purple",
    [1,1] => "orange"
);

function plotA(G::Vector{ComplexF64})
    ω = getOmegaRange() .+ 1.572
    A = -imag(G) / pi
    xr = [0, 5]
    yr = [0, 1]
    xl = "ω / t"
    yl = "A(ω)"
    ll = "J = " * string(info["coupling constant"]) * "t"
    fs = :box

    gr(size=(800,640),legend=false, markerstrokewidth=0)
    plot(ω, A,
        xlim = xr,
        ylim = yr,
        xlabel = xl,
        ylabel = yl,
        label = ll,
        framestyle = fs,
        alpha = 0.5,
        markersize = 3,
        left_margin = 5mm,
        bottom_margin = 5mm
    )
end

numstring = (x, n) -> map(x -> string(x .+ 10^n)[2:end], collect(x))

function plotSave(plt, name)
    png(plt, name)
end

function plotP(peaks)
    gr(size=(800,640),legend=false,markerstrokewidth=0)
    plt = scatter()
    origin = collect(keys(peaks))[1]
    for (m, data) in pairs(peaks)
        keyCollection = collect(keys(data))

        for key in keyCollection
            scatter!(plt,
                [key for _ in 1:length(data[key])],
                data[key] .- peaks[origin][key][1],
                xlim = [0.3,1],
                ylim = [0,2.5],
                xlabel = "J / t",
                ylabel = "ω / t",
                color = colors[m],
                markershape = :circle,
                markeralpha = 0.333,
                markersize = 3,
                framestyle = :box,
                legend = false
            )
        end
    end
    return plt
end

function plotPS(peaks)
    gr(size=(800,640),legend=false,markerstrokewidth=0)
    plots = []
    origin = collect(keys(peaks))[1]
    for (m, data) in pairs(peaks)
        keyCollection = collect(keys(data))

        plt = scatter()
        for key in keyCollection
            scatter!(plt,
                [key for _ in 1:length(data[key])],
                data[key] .- peaks[origin][key][1],
                xlim = [0.3,1],
                ylim = [0,2.5],
                xlabel = "J / t",
                ylabel = "ω / t",
                color = colors[m],
                markershape = :circle,
                markeralpha = 0.333,
                markersize = 3,
                framestyle = :box,
                legend = false
            )
        end
        push!(plots, plt)

    end
    return plots
end

function plotPS2(peaks)
    gr(size=(800,640),legend=false,markerstrokewidth=0)
    plots = []
    origin = collect(keys(peaks))[1]
    for (m, data) in pairs(peaks)
        keyCollection = collect(keys(data))

        plt = scatter()
        for key in keyCollection
            x = [key for _ in 1:length(data[key][1])]
            y = data[key][1] .- peaks[origin][key][1][1]
            w = data[key][2]

            scatter!(plt,
                x, y,
                xlim = [0.0,1],
                ylim = [0,5],
                xlabel = "J / t",
                ylabel = "ω / t",
                title = "[m₄, m₃] = " * string(m),
                color = colors[m],
                markershape = :circle,
                markeralpha = w,
                markersize = 2,
                framestyle = :box,
                legend = false,
                right_margin = 5mm,
                top_margin = 2mm,
                left_margin = 5mm,
                bottom_margin = 5mm,
                xtickfont = font(16),
                ytickfont = font(16),
                guidefont = font(16),
                titlefont = font(16)
            )
        end
        push!(plots, plt)

    end
    return plots
end


function plotS(peaks)
    gr(size=(800,640))
    plots = []
    origin = collect(keys(peaks))[1]

    for (m, data) in pairs(peaks)
        keyCollection = collect(keys(data))

        x = keyCollection
        ω = getOmegaRange()
        ωpts = length(ω)
        y = ω .- ω[1]
        z = zeros(Float64, length(y), length(x))
        for key in keyCollection
            xi = findfirst(q -> q == key, x)
            ωi = findfirst(q -> q == peaks[origin][key][1][1], ω)
            A = -imag(data[key][3]) / pi
            z[1:(ωpts - ωi + 1), xi] .= A[ωi:ωpts]
        end

        mLabels = [
            "\\mathrm{Bare}",
            "m_4 = " * string(m)[2:(end-1)],
            "m_4, m_3 = " * string(m)[2:(end-1)]
        ]

        xt = 0.0:0.2:1.0
        yt = 0:2:10
        fsize = 28
        plt = heatmap(
            x, y, z,
            xlim = [0.0,1],
            ylim = [0,10],
            xticks = (xt, [if isinteger(x) latexstring(Int64(x)) else latexstring(x) end for x in xt]),
            yticks = (yt, latexstring.(yt)),
            clim = (0,0.5),
            c = cgrad(:hot, rev=true),
            colorbar = true,
            xlabel = L"J / t",
            ylabel = latexstring("\\omega / t"),
            framestyle = :box,
            legend = false,
            title = latexstring(mLabels[1 + length(m)]),
            right_margin = 10mm,
            top_margin = 2mm,
            left_margin = 2mm,
            bottom_margin = 2mm,
            fontfamily = "times",
            xtickfontsize = fsize,
            ytickfontsize = fsize,
            guidefontsize = fsize,
            titlefontsize = fsize
        )
        push!(plots, plt)

    end
    return plots
end
