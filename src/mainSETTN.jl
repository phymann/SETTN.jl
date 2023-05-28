include("misc.jl")

function getlsTrHn(H::MPO, s; kwargs...)
    oplev = get(kwargs, :oplev, 2)
    ldTrHnQ = get(kwargs, :ldTrHnQ, false)
    nmax = get(kwargs, :nmax, 3)
    nrmHnQ = get(kwargs, :nrmHnQ, false)
    fname = get(kwargs, :fname, "test")

    fname *= "_TrHn.jld2"

    if ldTrHnQ
        file = jldopen(fname, "r")
        @unpack lsTrHn, nrm0 = file
    else
        H0 = MPO(s, "Id")
        nrm0 = tr(H0)
        H0 /= nrm0

        lsTrHn = []
        push!(lsTrHn, tr(H0))

        Hn = copy(H)
        Hn /= nrm0
        push!(lsTrHn, tr(Hn))

        nrmtot = 1.0
        for i = 2:nmax
            if oplev > 0
                println("cal tr(Hn) for n = $i, nmax = $nmax")
            end
            if nrmHnQ
                # this method may result in large numerical error!!!
                Hn, nrm = jwapply(H/i, Hn; kwargs...)
            else
                Hn, nrm = jwapply(H, Hn; kwargs...)
            end
            nrmtot *= nrm
            push!(lsTrHn, tr(Hn) * nrmtot)
        end
        file = jldopen(fname, "w")
        @pack! file = lsTrHn, nrm0
        @show fname
    end
    close(file)

    return lsTrHn, nrm0
end

function getFe(lsTrHn, β, nrm0; kwargs...)
    diftol = get(kwargs, :diftol, 1e-16)
    nrmHnQ = get(kwargs, :nrmHnQ, false)

    nmax = length(lsTrHn)
    fe0 = -β^-1 * log(nrm0)
    fe = fe0

    sumls = 0
    for i = 1:nmax
        if nrmHnQ
            sumls += lsTrHn[i]* (-β)^(i-1)
        else
            sumls += convert(Float64, lsTrHn[i]* (-big(β))^(i-1)/factorial(big(i-1)))
        end
        feold = fe
        fe = fe0 - β^-1 * log(sumls)
        if i>2 && abs((fe-feold)/feold) < diftol
            println("SETTN converges at n = $i")
            break
        end
        if i == nmax
            println("SETTN NOT converges even at n = $i")
        end
    end

    return fe
end

function getρ(H::MPO, s, β; kwargs...)
    oplev = get(kwargs, :oplev, 2)
    nmax = get(kwargs, :nmax, 3)
    nrmHnQ = get(kwargs, :nrmHnQ, false)

    H0 = MPO(s, "Id")
    nrm0 = tr(H0)
    H0 /= nrm0

    lsTrHn = []
    push!(lsTrHn, tr(H0))

    Hn = copy(H)
    Hn /= nrm0
    push!(lsTrHn, tr(Hn))

    nrmtot = 1.0
    for i = 2:nmax
        if oplev > 1
            println("cal tr(Hn) for n = $i / $nmax")
        end
        Hn, nrm = jwapply(H/i, Hn; kwargs...)
        nrmtot *= nrm
        push!(lsTrHn, tr(Hn) * nrmtot)
    end

    return lsTrHn, nrm0
end

function mainSETTN(H::MPO, s, lsβ; kwargs...)
    oplev = get(kwargs, :oplev, 0)
    lsTrHn, nrm0 = getlsTrHn(H, s; kwargs...)

    lsFE = copy(lsβ)
    nβ = length(lsβ)
    for (idx, β) in enumerate(lsβ)
        if oplev > 1
            println("cal Fe for β = $β, i = $idx / $nβ")
        end
        lsFE[idx] = getFe(lsTrHn, β, nrm0; kwargs...)
    end

    return lsFE
end
