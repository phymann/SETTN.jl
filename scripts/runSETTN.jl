using ITensors
using JLD2
using UnPack
using FlexiMaps
using LinearAlgebra
using MKL
using Infiltrator

include("src/mainSETTN.jl")
include("src/mpo.jl")
include("src/misc.jl")

function runSETTN(maxdim)
    println("maxdim = $maxdim")
    # ======================
    # setup multiple threads

    BLAS.set_num_threads(Sys.CPU_THREADS)
    ITensors.Strided.set_num_threads(1)
    ITensors.enable_threaded_blocksparse(false)

    # show info
    println()
    @show pwd()
    @show Threads.nthreads()
    @show Sys.CPU_THREADS
    @show BLAS.get_num_threads()
    @show ITensors.Strided.get_num_threads()
    @show ITensors.using_threaded_blocksparse()
    println()

    # ================
    # setup parameters
    nmax = 512
    cutoff = 1e-16
    SwpConvCrit = 1e-8
    symmQ = false
    L = 12
    hz = .42
    lsβ = LogRange(.1, 42, 8)
    si = 6
    fname = "rslt/SwpConvCrit=$(SwpConvCrit)_symm=$(symmQ)_maxdim$(maxdim)"
    if !isdir("rslt")
        mkdir("rslt")
    end
    opnames = ["Sz"]
    # save
    para = Dict{Symbol, Any}()
    @pack! para = nmax,
    cutoff,
    maxdim,
    symmQ,
    L,
    lsβ,
    fname,
    opnames,
    hz,
    si,
    SwpConvCrit

    # ==============
    # main functions
    s = siteinds("S=1/2", L; conserve_qns = symmQ)
    H = MPO(heisenberg(L, hz), s)

    lsfeED, lsexED = simpleED(L, hz, lsβ, si)

    # fe_settn
    lsfe, lsex = mainSETTN(H, s, lsβ, opnames;
        nmax = nmax,
        cutoff = cutoff,
        maxdim = maxdim,
        fname = fname,
        SwpConvCrit = SwpConvCrit
    )

    # ============
    # save results
    rslt = Dict{Symbol, Any}()
    @pack! rslt = lsex, lsexED, lsfe, lsfeED
    file = jldopen(fname*".jld2", "w")
    @pack! file = rslt, para
    close(file)
    @show fname
end
