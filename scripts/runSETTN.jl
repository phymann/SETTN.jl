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

    if read(`hostname`,String)[1:end-1] == "jw-mbp.local"
        MKLthreads = 1
        stridedThreads = Sys.CPU_THREADS
    else
        MKLthreads = Sys.CPU_THREADS
        stridedThreads = 1
    end
    BLAS.set_num_threads(MKLthreads)
    ITensors.Strided.set_num_threads(stridedThreads)
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
    symmQ = false
    ns = 8
    hz = .42
    lsβ = LogRange(.1, 42, 8)
    si = 6
    fname = "rslt/symm$(symmQ)_maxdim$(maxdim)"
    opnames = ["Sz"]
    # save
    para = Dict{Symbol, Any}()
    @pack! para = nmax,
    cutoff,
    maxdim,
    symmQ,
    ns,
    lsβ,
    fname,
    opnames,
    hz,
    si

    # ==============
    # main functions
    s = siteinds("S=1/2", ns; conserve_qns = symmQ)
    H = MPO(heisenberg(ns, hz), s)

    lsfeED, lsexED = simpleED(ns, hz, lsβ, si)

    # fe_settn
    lsfe, lsex = mainSETTN(H, s, lsβ, opnames;
        nmax = nmax,
        cutoff = cutoff,
        maxdim = maxdim,
        fname = fname
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
