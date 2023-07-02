using SETTN.ITensors

include("misc.jl")

@testset "SETTN for Heisenberg chian" begin
    lsβ = LogRange(1,64,4)
    ns = 6
    s = siteinds("S=1/2", ns; conserve_qns = true)
    os = OpSum()
    for j = 1:(ns-1)
        os += 0.5, "S+", j, "S-", j+1
        os += 0.5, "S-", j, "S+", j+1
        os += 1.0, "Sz", j, "Sz", j+1
    end
    H = MPO(os, s)

    lsfe = zeros(size(lsβ))
    # lsex = zeros(size(lsβ))
    for (βi, β) in enumerate(lsβ)
        # density matrix at β/2: ρ(β/2)
        rho1, nrm0 = getρ(H, β/2;
            maxord = 1024,
            cutoff = 1e-16,
            maxdim = 512,
            outputlevel = 0,
            SwpConvCrit = 1e-24
        )
        lsfe[βi] = getFe(rho1, β, nrm0)
    end

    lsfe_ed = ED_Heisenberg(ns, lsβ)

    @test maximum(abs.((lsfe - lsfe_ed)./lsfe_ed)) < 1e-6
end
