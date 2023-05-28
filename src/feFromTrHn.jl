using ITensors
using JLD2
using UnPack
using FlexiMaps
using CairoMakie
using MakiePublication

include("mainSETTN.jl")
include("mainED.jl")

# setup multiple threads
# ----------------------
MKLthreads = 1
BLAS.set_num_threads(MKLthreads)
ITensors.Strided.set_num_threads(Sys.CPU_THREADS)
ITensors.enable_threaded_blocksparse(false)

println("-----------------------------------------------------------")
@show pwd()
@show Threads.nthreads()
@show Sys.CPU_THREADS
@show BLAS.get_num_threads()
@show ITensors.Strided.get_num_threads()
@show ITensors.using_threaded_blocksparse()
println()

ns = 12
s = siteinds("S=1/2", 12; conserve_qns = true)

nmax = 128

function heisenberg(n)
    os = OpSum()
    for j = 1:(n-1)
        os += 0.5, "S+", j, "S-", j+1
        os += 0.5, "S-", j, "S+", j+1
        os += 1.0, "Sz", j, "Sz", j+1
    end
    return os
end

H = MPO(heisenberg(ns), s)

lsβ = collect(maprange(log, 1e-3, 1e-1, length=32))
nβ = length(lsβ)

fe_ed, _ = mainED(H, s, lsβ; loadQ = true)

fname = "rslt/nmax129TrHn.jld2"
file = jldopen(fname, "r")
@unpack lsTrHn, nrm0 = file
close(file)

fe_settn = copy(lsβ)
for (idx, β) in enumerate(lsβ)
    println("cal Fe for βi = $idx / $nβ")
    fe_settn[idx] = getFe(lsTrHn, β, nrm0;
    tol_fediff = 1e-24)
end

fname = "rslt/nmax$nmax.jld2"
file = jldopen(fname, "w")
@pack! file = fe_settn, fe_ed, lsβ
close(file)
@show fname
