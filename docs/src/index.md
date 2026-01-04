# QSymbolic.jl

*A Julia package for symbolic quantum mechanics with basis-aware state representations.*

## Features

- **Hilbert Spaces**: Define finite and infinite-dimensional Hilbert spaces
- **Explicit Bases**: Named orthonormal bases with automatic orthogonality  
- **Flexible Indices**: Symbolic, numeric, or multi-index kets for composite systems
- **Basis Transforms**: Register transformations between bases; cross-basis inner products computed automatically
- **Composite Systems**: Tensor products with order-independent (bosonic) behavior and factorized transforms
- **Custom Contraction Rules**: Define non-orthonormal inner products for dressed states
- **Operators**: Outer product operators `|ψ⟩⟨ϕ|`, operator algebra, and function-defined operators
- **Symbolics.jl Backend**: Full symbolic computation with `Sym`, `@variables`, `KroneckerDelta`

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/TendonFFF/QSymbolic.jl")
```

## Quick Start

```julia
using QSymbolic

# Create a 2-dimensional Hilbert space with default basis
H, Hb = HilbertSpace(:spin, 2)

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

# Symbolic variables (Symbolics.jl backend)
n = Sym(:n, :nonnegative, :integer)
√n + 1  # → 1 + √n

# Symbolic inner products
F, Fb = FockSpace(:mode)
ket_n = Ket(Fb, n)
m = Sym(:m, :nonnegative, :integer)
Ket(Fb, m)' * ket_n  # → δ(m,n) (KroneckerDelta)
```

## Type Hierarchy

### States
```
AbstractKet
├── Ket           # Basic ket with index: |ψ⟩
├── ProductKet    # Tensor product: |ψ⟩⊗|ϕ⟩ (order-independent)
├── WeightedKet   # Scalar × ket: α|ψ⟩
└── SumKet        # Superposition: α|ψ⟩ + β|ϕ⟩

AbstractBra       # Lazy adjoints of kets
├── Bra, ProductBra, WeightedBra, SumBra
```

### Operators
```
AbstractOperator
├── Outer            # Single |ψ⟩⟨ϕ|
├── Operator         # Sum of weighted outers
├── Identity         # Identity 𝕀
└── FunctionOperator # User-defined action
```

## Contents

```@contents
Pages = [
    "guide/getting_started.md",
    "guide/transforms.md",
    "guide/composite.md",
    "guide/operators.md",
    "guide/symbolic.md",
    "guide/contraction_rules.md",
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
