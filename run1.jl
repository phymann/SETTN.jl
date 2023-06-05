include("src/runSETTN.jl")
@time runSETTN.([2^i for i in 6:9])
