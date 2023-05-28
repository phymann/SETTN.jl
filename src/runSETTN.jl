using ITensors
using JLD2
using UnPack
using FlexiMaps
using LinearAlgebra
using MKL

include("mainSETTN.jl")
include("mainED.jl")
include("ham.jl")
include("misc.jl")

let 
    for maxdim = Int64.(floor.(map(x->2^x, [6:8;9.5])))
        println("maxdim = $maxdim")
    # setup multiple threads
    # --------------------------------------------------------------------------------------
    MKLthreads = 8
    if MKLthreads == 1
        stridedThreads = Sys.CPU_THREADS
    else
        stridedThreads = 1
    end
    BLAS.set_num_threads(MKLthreads)
    ITensors.Strided.set_num_threads(stridedThreads)
    ITensors.enable_threaded_blocksparse(false)

    println()
    @show pwd()
    @show Threads.nthreads()
    @show Sys.CPU_THREADS
    @show BLAS.get_num_threads()
    @show ITensors.Strided.get_num_threads()
    @show ITensors.using_threaded_blocksparse()
    println()
    # --------------------------------------------------------------------------------------

    # setup parameters
    # --------------------------------------------------------------------------------------
    ldTrHnQ = false
    ldEDQ = true
    nmax = 256
    nrmHnQ = false
    cutoff = 1e-16
    toldif = 1e-16
    # maxdim = maxdim
    @show maxdim
    oplev = 42
    symmQ = false
    ns = 12
    lsβ = LogRange(0.42, 42, 32)
    fname = "rslt/nmax$(nmax)_symm$(symmQ)_maxdim$(maxdim)"
    # save them all
    para = Dict{Symbol, Any}()
    @pack! para = nmax, nrmHnQ, cutoff, toldif, maxdim, oplev, symmQ, ns, lsβ, fname, ldEDQ
    # para = NamedTuple(para)
    # --------------------------------------------------------------------------------------

    # main functions
    # --------------------------------------------------------------------------------------
    s = siteinds("S=1/2", ns; conserve_qns = symmQ)
    H = MPO(heisenberg(ns), s)

    fe_ed, _ = mainED(H, s, lsβ;
        ldEDQ = ldEDQ,
        # fname = fname
    )

    fe_settn = mainSETTN(H, s, lsβ;
        nmax = nmax,
        oplev = oplev,
        nrmHnQ = nrmHnQ,
        ldTrHnQ = ldTrHnQ,
        cutoff = cutoff,
        toldif = toldif,
        maxdim = maxdim,
        fname = fname
    )
    # --------------------------------------------------------------------------------------

    # save results
    # --------------------------------------------------------------------------------------
    rslt = Dict{Symbol, Any}()
    @pack! rslt = fe_settn, fe_ed
    # rslt = NamedTuple(rslt)

    file = jldopen(fname*".jld2", "w")
    @pack! file = rslt, para
    close(file)
    @show fname
    # --------------------------------------------------------------------------------------
end
end
