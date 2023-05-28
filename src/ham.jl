function heisenberg(n)
    os = OpSum()
    for j = 1:(n-1)
        os += 0.5, "S+", j, "S-", j+1
        os += 0.5, "S-", j, "S+", j+1
        os += 1.0, "Sz", j, "Sz", j+1
    end
    return os
end
