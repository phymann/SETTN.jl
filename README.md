# SETTN.jl

The series-expansion thermal tensor network (SETTN) algorithm [[PRB **95**, 161104](https://doi.org/10.1103/PhysRevB.95.161104)], for quantum many-body lattice models, is implemented.

## Installation

```julia
] add SETTN
```

## Usage

There are only two functions exported,

```julia
"""
    input: β and H
    output: ρ(β) and norm0
"""
function getρ(H::MPO, β::Float64; kwargs...)
```

and

```julia
"""
    input: ρ(β/2), β and norm0 from `getρ`
    output: Fe(β)
"""
function getFe(rho::MPO, β::Float64, nrm0::Float64)
```


One may check the example given in the `test` folder.

> Note This package is mainly used for other more advanced thermal tensor network algorithm, such as XTRG (to be implemented) and [tanTRG.jl](https://github.com/phymann/tanTRG.jl).

## Details

Two-site variational MPO sum and product, such as the one used in [PRB **95**, 161104]((https://doi.org/10.1103/PhysRevB.95.161104)) and [PRX **8**, 031082](https://doi.org/10.1103/PhysRevX.8.031082), is implemented.

## TODO

- Implement a more numerically stable MPO product, such as the one given in [PRB **102**, 035147](https://doi.org/10.1103/PhysRevB.102.035147)
- Make two-site variational MPO sum and product available for [`ITensors.jl`](https://github.com/ITensor/ITensors.jl) package
