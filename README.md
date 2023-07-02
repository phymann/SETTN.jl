# SETTN.jl

This repo implements [Series-expansion thermal tensor network (SETTN)](https://doi.org/10.1103/PhysRevB.95.161104) approach for quantum many-body lattice models.

## How to use

After cloning this repo, from its main directory, start julia using

```shell
julia -t 8 --project=.
```

Then

```julia
using Pkg; Pkg.instantiate()
```

Finally, start from `src/runSETTN.jl`, where spin-1/2 Heisenberg chain is calculated.

Other models can be studied by constructing your own Hamiltonian, via [`ITensors.jl`](https://github.com/ITensor/ITensors.jl) package.

## TODO

- [DONE] ~~Implement two-site variational MPO product, such as the one used in [PRB **95**, 161104]((https://doi.org/10.1103/PhysRevB.95.161104)) and [PRX **8**, 031082](https://doi.org/10.1103/PhysRevX.8.031082)~~
- Implement a more numerically stable MPO product, such as the one given in [PRB **102**, 035147](https://doi.org/10.1103/PhysRevB.102.035147)
- Make this repo a standard Julia package
