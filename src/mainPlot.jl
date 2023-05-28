using CairoMakie
using MakiePublication
using JLD2
using UnPack

struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

let
    # --------------------------------------------------------------------------------------
    symmQ = true
    nmax = 256
    maxdim = 512
    fname = "rslt/NEWnmax$(nmax)_symm$(symmQ)_maxdim$maxdim.jld2"

    file = jldopen(fname, "r")
    @unpack fe_settn, fe_ed, lsβ = file
    close(file)
    x = lsβ.^-1
    y = abs.((fe_settn-fe_ed)./fe_ed)
    y .+= 1e-16

    # --------------------------------------------------------------------------------------
    symmQ = false
    nmax = 256
    maxdim = 512
    fname = "rslt/NEWnmax$(nmax)_symm$(symmQ)_maxdim$maxdim.jld2"

    file = jldopen(fname, "r")
    @unpack fe_settn, fe_ed, lsβ = file
    close(file)
    x1 = lsβ.^-1
    y1 = abs.((fe_settn-fe_ed)./fe_ed)
    y1 .+= 1e-16

    # --------------------------------------------------------------------------------------
    symmQ = true
    nmax = 128
    maxdim = 512
    fname = "rslt/NEWnmax$(nmax)_symm$(symmQ)_maxdim$maxdim.jld2"

    file = jldopen(fname, "r")
    @unpack fe_settn, fe_ed, lsβ = file
    close(file)
    x2 = lsβ.^-1
    y2 = abs.((fe_settn-fe_ed)./fe_ed)
    y2 .+= 1e-16

    # --------------------------------------------------------------------------------------
    symmQ = true
    nmax = 256
    maxdim = 1024
    fname = "rslt/NEWnmax$(nmax)_symm$(symmQ)_maxdim$maxdim.jld2"

    file = jldopen(fname, "r")
    @unpack fe_settn, fe_ed, lsβ = file
    close(file)
    x2 = lsβ.^-1
    y2 = abs.((fe_settn-fe_ed)./fe_ed)
    y2 .+= 1e-16

    # # --------------------------------------------------------------------------------------
    # symmQ = false
    # nmax = 256
    # maxdim = 1024
    # fname = "rslt/nmax$(nmax)_symm$(symmQ)_maxdim$maxdim.jld2"

    # file = jldopen(fname, "r")
    # @unpack fe_settn, fe_ed, lsβ = file
    # close(file)
    # x3 = lsβ.^-1
    # y3 = abs.((fe_settn-fe_ed)./fe_ed)
    # y3 .+= 1e-16

    function plot1()
        f = Figure(
            fontsize = 24,
            resolution = (600, 600)
        )

        ax = Axis(
            f[1, 1];
            title = "1D Heisenberg model",
            yscale = log10,
            ylabel = L"(F_{\text{SETTN}} - F_{\text{ED}})/F_{\text{ED}}",
            xscale = log10,
            xlabel = L"T",
            xminorticks = IntervalsBetween(10),
            yminorticks = IntervalsBetween(10),
            xticks = LogTicks(IntegerTicks()),
            # yticks = LogTicks(IntegerTicks())
        )

        sl1 = scatterlines!(ax, x, y;
            marker = :vline,
            markersize = 24
        )
        sl2 = scatterlines!(ax, x1, y1;
            marker = :circle,
            markersize = 8
        )
        sl3 = scatterlines!(ax, x2, y2;
            marker = :hline,
            markersize = 24
        )
        sl4 = scatterlines!(ax, x2, y2;
            marker = :diamond,
            markersize = 12
        )
        # sl4 = scatterlines!(ax, x3, y3; )
        axislegend(
            ax,
            [
                sl1,
                sl2,
                sl3,
                sl4
            ],
            [
                L"D = 512,\quad n_{\text{max}} = 256,\quad \text{ symm}",
                L"D = 512,\quad n_{\text{max}} = 256",
                L"D = 512,\quad n_{\text{max}} = 128,\quad \text{ symm}",
                L"D = 1024,\quad n_{\text{max}} = 128,\quad \text{ symm}",
            ]
        )

        ylims!(ax, 1e-16, 100)
        current_figure()
    end
    with_theme(plot1, theme_web())
end
