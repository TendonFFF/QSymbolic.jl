# QSymbolic.jl

*A Julia package for symbolic quantum mechanics with basis-aware state representations.*

## Features

- **Hilbert Spaces**: Define finite and infinite-dimensional Hilbert spaces
- **Explicit Bases**: Named orthonormal bases with automatic orthogonality
- **Basis Transforms**: Register transformations between bases; cross-basis inner products computed automatically
- **Composite Systems**: Tensor products with factorized transforms
- **Operators**: Outer product operators `|ψ⟩⟨ϕ|`, operator algebra, and function-defined operators
- **Symbolic Scalars**: Lazy arithmetic with symbolic variables for truly symbolic computation
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
Hb = Basis(H, :default)

# Define spin-z basis
Zb = Basis(H, :z)

# Create kets in the z-basis
up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Orthonormality
up' * up    # → 1
up' * down  # → 0

# Build operators via outer products
P_up = up * up'      # |↑⟩⟨↑| (projector)
σ_plus = up * down'  # |↑⟩⟨↓| (raising operator)

# Apply operators
P_up * up    # → |↑⟩
P_up * down  # → 0

# Pauli Z operator
σz = P_up - (down * down')
σz * up    # → |↑⟩
σz * down  # → -|↓⟩

# Symbolic scalars for lazy evaluation
n = Sym(:n)
expr = √n * 2
substitute(expr, :n => 4) |> evaluate  # → 4.0
```

## Contents

```@contents
Pages = [
    "guide/getting_started.md",
    "guide/transforms.md",
    "guide/composite.md",
    "guide/operators.md",
    "guide/symbolic.md",
    "api/spaces.md",
    "api/bases.md",
    "api/states.md",
    "api/transforms.md",
    "api/operators.md",
    "api/symbolic.md",
]
Depth = 2
```

## Index

```@index
```
