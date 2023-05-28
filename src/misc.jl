macro nt(ex)
    Expr(:tuple, [Expr(:(=), esc(arg), arg) for arg in ex.args]...)
end

Base.NamedTuple(d::Dict{Symbol, Any}) = NamedTuple{Tuple(keys(d))}(values(d))

"""
    normalize the local tensor after each svd, also return the total normalization factor
"""
function jwcontract(A::MPO, B::MPO; kwargs...)
    if hassameinds(siteinds, A, B)
        error(
            "In `contract(A::MPO, B::MPO)`, MPOs A and B have the same site indices. The indices of the MPOs in the contraction are taken literally, and therefore they should only share one site index per site so the contraction results in an MPO. You may want to use `replaceprime(contract(A', B), 2 => 1)` or `apply(A, B)` which automatically adjusts the prime levels assuming the input MPOs have pairs of primed and unprimed indices.",
        )
    end
    cutoff::Float64 = get(kwargs, :cutoff, 1e-14)
    resp_degen::Bool = get(kwargs, :respect_degenerate, true)
    maxdim::Int = get(kwargs, :maxdim, maxlinkdim(A) * maxlinkdim(B))
    mindim::Int = max(get(kwargs, :mindim, 1), 1)
    N = length(A)
    N != length(B) &&
        throw(DimensionMismatch("lengths of MPOs A ($N) and B ($(length(B))) do not match"))
    # Special case for a single site
    N == 1 && return MPO([A[1] * B[1]])
    A = ITensors.orthogonalize(A, 1)
    B = ITensors.orthogonalize(B, 1)
    A = sim(linkinds, A)
    sA = siteinds(uniqueinds, A, B)
    sB = siteinds(uniqueinds, B, A)
    C = MPO(N)
    lCᵢ = Index[]
    R = ITensor(true)
    nrm = 1.0
    for i in 1:(N - 2)
        RABᵢ = R * A[i] * B[i]
        left_inds = [sA[i]..., sB[i]..., lCᵢ...]
        C[i], S, R = ITensors.svd(
            RABᵢ,
            left_inds;
            # ortho="left",
            # tags=commontags(linkinds(A, i)),
            lefttags=ITensors.commontags(linkinds(A, i)),
            cutoff=cutoff,
            maxdim=maxdim,
            mindim=mindim,
            kwargs...,
        )
        nrmS = norm(matrix(S))
        nrm *= nrmS
        S /= nrmS
        R *= S
        lCᵢ = dag(commoninds(C[i], R))
    end
    # ------------------------------------------------------
    # special care should be taken for the penultimate site!
    i = N - 1
    RABᵢ = R * A[i] * B[i] * A[i + 1] * B[i + 1]
    left_inds = [sA[i]..., sB[i]..., lCᵢ...]
    C[N - 1], S, C[N] = ITensors.svd(
    RABᵢ,
    left_inds;
    # ortho="right",
    # tags=commontags(linkinds(A, i)),
    righttags=ITensors.commontags(linkinds(A, i)),
    cutoff=cutoff,
    maxdim=maxdim,
    mindim=mindim,
    kwargs...,
    )
    nrmS = norm(matrix(S))
    nrm *= nrmS
    S /= nrmS
    C[N-1] *= S
    truncate!(C; kwargs...)
    return C, nrm
end


function jwapply(A::MPO, B::MPO; kwargs...)
    AB, nrm = jwcontract(A', B; kwargs...)
    return replaceprime(AB, 2 => 1), nrm
end
