using LinearAlgebra
using FlexiMaps

LogRange(ls1,lsend,nls) = maprange(log, ls1, lsend, length=nls)

"""
    return free energies at a given list of β points
"""
function ED_Heisenberg(ns, lsβ)
    mz = [1/2 0; 0 -1/2]
    mx = [0 1/2; 1/2 0]
    # my = [0 -1im/2; 1im/2 0]
    my1 = [0 1/2; -1/2 0]
    m1 = [1 0; 0 1]

    # create full matrices
    full1(m,i) = kron([idx == i ? m : m1 for idx in 1:ns]...)
    full2(ma,mb,ia,ib) = kron([idx == ia ? ma : (idx == ib ? mb : m1) for idx in 1:ns]...)

    ham = zeros(Float64, 2^ns,2^ns)
    for i = 1:ns-1
        ham += full2(mx,mx,i,i+1)
        ham += -full2(my1,my1,i,i+1)
        ham += full2(mz,mz,i,i+1)
    end
    # for i = 1:ns
    #     ham += -hz * full1(mz,i)
    # end
    vals, _ = eigen(ham)

    # sz6 = full1(mz,si)
    # lsex = copy(lsβ)
    lsfe = copy(lsβ)

    for (βi,β) in enumerate(lsβ)
        bigz = tr(exp(-β*ham))
        bigz = sum(exp.(-β .* vals))
        lsfe[βi] = -1/β * log(bigz)
        # lsex[βi] = sum([vecs[:,i]' * sz6 * vecs[:,i] / (vecs[:,i]' * vecs[:,i]) *
        #  exp(-β*vals[i]) for i in eachindex(vals)]) / bigz
    end

    return lsfe
end

function ED_ising(ns, hx, lsβ, si)
    mz = [1/2 0; 0 -1/2]
    mx = [0 1/2; 1/2 0]
    m1 = [1 0; 0 1]

    # create full matrices
    full1(m,i) = kron([idx == i ? m : m1 for idx in 1:ns]...)
    full2(ma,mb,ia,ib) = kron([idx == ia ? ma : (idx == ib ? mb : m1) for idx in 1:ns]...)

    ham = zeros(Float64, 2^ns,2^ns)
    for i = 1:ns-1
        ham += -full2(mz,mz,i,i+1)
    end
    for i = 1:ns
        ham += -hx * full1(mx,i)
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

    return lsfe, lsex
end
