"""
    input: β and H
    output: ρ(β)
"""
function getρ(H::MPO, β::Float64; kwargs...)
    maxord = get(kwargs, :maxord, 1024)
    exactmax = get(kwargs, :exactmax, 3)
    stop_tol = get(kwargs, :stop_tol, 1e-16)
    s = firstsiteinds(H)
    s = dag(s) # for symm case

    H0 = MPO(s, "Id")
    nrm0 = norm(H0)
    H0 /= nrm0
    rho1 = H0

    Hn = copy(H)
    Hn /= nrm0

    rho1 = +(rho1, -β * Hn; alg="directsum")

    lognrmtot = 0
    for i = 2:maxord
        if i ≤ exactmax
            Hn = apply(H, Hn; alg="zipup", kwargs...)
        else
            Hn = apply(H, Hn; alg="variational", kwargs...)
        end
        nrm = norm(Hn)
        Hn /= nrm
        lognrmtot += log(nrm)
        coeff = Float64(((-1)^i*exp(lognrmtot + (i)*log(big(β)) - log(factorial(big(i))))))
        # @show coeff
        nrm1 = norm(rho1)
        nrm2 = coeff
        if i ≤ exactmax
            rho1 = +(rho1, coeff*Hn; alg="directsum")
        else
            # NB! The second MPO is expected to be closer to the resulting MPO
            if nrm1 > nrm2
                rho1 = +(coeff/nrm1*Hn, rho1/nrm1; alg="variational", kwargs...)
                rho1 *= nrm1
            else
                rho1 = +(rho1/nrm2, Hn; alg="variational", kwargs...)
                rho1 *= nrm2
            end
        end
        nrmnew = norm(rho1)
        # @show nrmnew
        stopQ = abs((nrmnew-nrm1)/nrm1)
        if stopQ < stop_tol
            println("ρ converged at $i/$maxord")
            break
        end
        if i == maxord
            println("ρ NOT converged")
        end
    end
    return rho1, nrm0
end

"""
    input: ρ(β/2) and β
    output: Fe(β)
"""
function getFe(rho::MPO, β::Float64, nrm0::Float64)
    return -1/β * (2*log(nrm0) + 2*log(norm(rho)))
end
