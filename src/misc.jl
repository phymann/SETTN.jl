using ITensors

macro nt(ex)
    Expr(:tuple, [Expr(:(=), esc(arg), arg) for arg in ex.args]...)
end

Base.NamedTuple(d::Dict{Symbol, Any}) = NamedTuple{Tuple(keys(d))}(values(d))

# """
#     normalize the local tensor after each svd, also return the total normalization factor
# """
# function jwcontract(A::MPO, B::MPO; kwargs...)
#     if hassameinds(siteinds, A, B)
#         error(
#             "In `contract(A::MPO, B::MPO)`, MPOs A and B have the same site indices. The indices of the MPOs in the contraction are taken literally, and therefore they should only share one site index per site so the contraction results in an MPO. You may want to use `replaceprime(contract(A', B), 2 => 1)` or `apply(A, B)` which automatically adjusts the prime levels assuming the input MPOs have pairs of primed and unprimed indices.",
#         )
#     end
#     cutoff::Float64 = get(kwargs, :cutoff, 1e-14)
#     maxdim::Int = get(kwargs, :maxdim, maxlinkdim(A) * maxlinkdim(B))
#     mindim::Int = max(get(kwargs, :mindim, 1), 1)
#     N = length(A)
#     N != length(B) &&
#         throw(DimensionMismatch("lengths of MPOs A ($N) and B ($(length(B))) do not match"))
#     # Special case for a single site
#     N == 1 && return MPO([A[1] * B[1]])
#     A = ITensors.orthogonalize(A, 1)
#     B = ITensors.orthogonalize(B, 1)
#     A = sim(linkinds, A)
#     sA = siteinds(uniqueinds, A, B)
#     sB = siteinds(uniqueinds, B, A)
#     C = MPO(N)
#     lCᵢ = Index[]
#     R = ITensor(true)
#     nrm = 1.0
#     for i in 1:(N - 2)
#         RABᵢ = R * A[i] * B[i]
#         left_inds = [sA[i]..., sB[i]..., lCᵢ...]
#         C[i], S, R = ITensors.svd(
#             RABᵢ,
#             left_inds;
#             lefttags=ITensors.commontags(linkinds(A, i)),
#             cutoff=cutoff,
#             maxdim=maxdim,
#             mindim=mindim,
#             alg="recursive",
#         )
#         nrmS = norm(matrix(S))
#         nrm *= nrmS
#         S /= nrmS
#         R *= S
#         lCᵢ = dag(commoninds(C[i], R))
#     end
#     # ------------------------------------------------------
#     # special care should be taken for the penultimate site!
#     i = N - 1
#     RABᵢ = R * A[i] * B[i] * A[i + 1] * B[i + 1]
#     left_inds = [sA[i]..., sB[i]..., lCᵢ...]
#     C[N - 1], S, C[N] = ITensors.svd(
#     RABᵢ,
#     left_inds;
#     righttags=ITensors.commontags(linkinds(A, i)),
#     cutoff=cutoff,
#     maxdim=maxdim,
#     mindim=mindim,
#     alg="recursive"
#     )
#     nrmS = norm(matrix(S))
#     nrm *= nrmS
#     S /= nrmS
#     C[N-1] *= S
#     truncate!(C; kwargs...)
#     return C, nrm
# end

# function jwapply(A::MPO, B::MPO; kwargs...)
#     AB, nrm = jwcontract(A', B; kwargs...)
#     return replaceprime(AB, 2 => 1), nrm
# end

function heisenberg(n::Int64, hz::Real)
    os = OpSum()
    for j = 1:(n-1)
        os += 0.5, "S+", j, "S-", j+1
        os += 0.5, "S-", j, "S+", j+1
        os += 1.0, "Sz", j, "Sz", j+1
    end
    if !(hz ≈ 0)
        for j = 1:n
            os += -hz, "Sz", j
        end
    end
    return os
end

using LinearAlgebra
using JLD2

function simpleED(ns, hz, lsβ, si)
    # local op
    mz = [1/2 0; 0 -1/2]
    mx = [0 1/2; 1/2 0]
    my = [0 -1im/2; 1im/2 0]
    m1 = [1 0; 0 1]

    # create full matrices
    full1(m,i) = kron([idx == i ? m : m1 for idx in 1:ns]...)
    full2(ma,mb,ia,ib) = kron([idx == ia ? ma : (idx == ib ? mb : m1) for idx in 1:ns]...)

    ham = zeros(ComplexF64, 2^ns,2^ns)
    for i = 1:ns-1
        ham += full2(mx,mx,i,i+1)
        ham += full2(my,my,i,i+1)
        ham += full2(mz,mz,i,i+1)
    end
    for i = 1:ns
        ham += -hz * full1(mz,i)
    end
    vals, vecs = eigen(ham)

    sz6 = full1(mz,si)
    lsex = copy(lsβ)
    lsfe = copy(lsβ)

    for (βi,β) in enumerate(lsβ)
        bigz = tr(exp(-β*ham))
        bigz = sum(exp.(-β .* vals))
        lsfe[βi] = -1/β * log(bigz)
        lsex[βi] = sum([vecs[:,i]' * sz6 * vecs[:,i] / (vecs[:,i]' * vecs[:,i]) *
         exp(-β*vals[i]) for i in eachindex(vals)]) / bigz
    end

    jldsave("rslt/ED_L=$(ns)_hz=$(hz)_si=$(si).jld2"; lsfeED=lsfe,  lsexED=lsex)

    return lsfe, lsex
end
