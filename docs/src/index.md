# QSymbolic.jl

*A Julia package for symbolic quantum mechanics with basis-aware state representations.*

## Features

- **Hilbert Spaces**: Define finite and infinite-dimensional Hilbert spaces
- **Explicit Bases**: Named orthonormal bases with automatic orthogonality
- **Basis Transforms**: Register transformations between bases; cross-basis inner products computed automatically
- **Composite Systems**: Tensor products with factorized transforms
- **Symbolic Computation**: Unevaluated inner products when transforms are undefined

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/TendonFFF/QSymbolic.jl")
```

## Quick Start

```julia
using QSymbolic

# Create a 2-dimensional Hilbert space
H = HilbertSpace(:spin, 2)

# Define spin-z and spin-x bases
Zb = Basis(H, :z)
Xb = Basis(H, :x)

# Create kets in the z-basis
up = BasisKet(Zb, :↑)
down = BasisKet(Zb, :↓)

# Orthonormality in same basis
up' * up    # → 1
up' * down  # → 0

# Define how x-basis relates to z-basis
define_transform!(Xb, Zb) do idx
    idx == :↑ ? (up + down) / √2 : (up - down) / √2
end

# Cross-basis inner products now work
up_x = BasisKet(Xb, :↑)
up' * up_x  # → 1/√2
```

## Contents

```@contents
Pages = [
    "guide/getting_started.md",
    "guide/transforms.md",
    "guide/composite.md",
    "api/spaces.md",
    "api/bases.md",
    "api/states.md",
    "api/transforms.md",
]
Depth = 2
```

## Index

```@index
```
