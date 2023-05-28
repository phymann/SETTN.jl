using ITensors
using KrylovKit
using LinearAlgebra
using JLD2
using UnPack
# using MKL

include("fuse_inds.jl")

function mainED(H, s, lsβ; fuse=true, binary=true, ldEDQ = false, fname = "rslt/testED")
    nsite = length(s)
    if nsite > 16
        @warn "System size of $nsite is likely too large for exact diagonalization."
    end

    fname *= "_vals.jld2"

    if ldEDQ
        file = jldopen(fname, "r")
        @unpack vals = file
    else
        if fuse
            if binary
            println("Fuse the indices using a binary tree")
            T = fusion_tree_binary(s)
            H_full = @time fuse_inds_binary(H, T)
            else
            println("Fuse the indices using an unbalances tree")
            T = fusion_tree(s)
            H_full = @time fuse_inds(H, T)
            end
        else
            println("Don't fuse the indices")
            @disable_warn_order begin
            H_full = @time contract(H)
            end
        end
        vals = eigvals(array(H_full))

        file = jldopen(fname, "w")
        @pack! file = vals
        @show fname
    end
    close(file)

    lsFe = copy(lsβ)
    lsIe = copy(lsβ)
    for (cnt, β) in enumerate(lsβ)
        bigz = sum(exp.(-β*vals))
        lsFe[cnt] = -β^-1 * log(bigz)
        lsIe[cnt] = sum(vals .* exp.(-β*vals))/bigz
    end

    return lsFe, lsIe
end
