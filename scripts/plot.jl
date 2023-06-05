using JLD2
using CairoMakie
using MakiePublication
using UnPack
using Infiltrator

function jwf(lsmaxdim)
    si = 6
    function jwp()
        f = Figure(
            resolution = (1200, 600),
            fontsize = 24
        )
        ax = Axis(f[1,1],
        ylabel = L"\langle S^z(i=%$(si)) \rangle",
        yscale = log10,
        xscale = log10
        )
        ax1 = Axis(f[1,2],
        yscale = log10,
        xscale = log10,
        ylabel = L"\text{relative err of free energy}"
        )
        for maxdim in lsmaxdim
            fn = "rslt/L=12_S=false_MD=$maxdim"
            @show fn
            file = jldopen(fn*".jld2","r")
            @unpack lsex, lsexED, lsfe, lsfeED = file["rslt"]
            @unpack  lsβ = file["para"]
            close(file)
            lsex = map(x->x[1][si], lsex)

            scatterlines!(ax, lsβ.^-1, abs.((lsex-lsexED)./lsexED),
            label = L"D=%$maxdim"
            )


            scatterlines!(ax1, lsβ.^-1, abs.((lsfe-lsfeED)./lsfeED),
            label = L"D=%$maxdim"
            )
        end
        axislegend(ax, position = :rt)
        axislegend(ax1, position = :lb)

        current_figure()
    end
    fig = with_theme(jwp, theme_web())
    # savefig("plots/L12.pdf", fig)
end
