using ITensors

macro nt(ex)
    Expr(:tuple, [Expr(:(=), esc(arg), arg) for arg in ex.args]...)
end

Base.NamedTuple(d::Dict{Symbol, Any}) = NamedTuple{Tuple(keys(d))}(values(d))

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
