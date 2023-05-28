lsbdim = Int64.(floor.(map(x->2^x, [6:8;9.5])))
for maxdim = lsbdim
    println("maxdim = $maxdim")
    include("runSETTN.jl")
end
