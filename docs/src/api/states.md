# States

Quantum states are represented as kets `|ψ⟩` and their duals as bras `⟨ψ|`. QSymbolic.jl provides types for basis states, superpositions, and product states.

## Overview

| Type | Description |
|:-----|:------------|
| `Ket` / `Bra` | Single basis state |
| `WeightedKet` / `WeightedBra` | Scalar × basis state |
| `SumKet` / `SumBra` | Linear combination of basis states |
| `ProductKet` / `ProductBra` | Tensor product of states (order-independent) |

## Index Types

Kets support flexible indices:

```julia
# Single index
Ket(basis, :ψ)           # Symbol
Ket(basis, 0)            # Integer (Fock states)
Ket(basis, Sym(:n))      # Symbolic

# Multi-index (for composite bases)
Ket(basis, (n, m))       # Tuple of indices
Ket(basis, (n, m, k))    # Arbitrary length
```

## Single-System States

### Abstract Types

```@docs
AbstractKet
AbstractBra
```

### Basis States

```@docs
Ket
Bra
```

### Weighted States

```@docs
WeightedKet
WeightedBra
```

### Superpositions

```@docs
SumKet
SumBra
```

## Composite States

For multi-particle systems. `ProductKet` is **order-independent** (bosonic).

```@docs
ProductKet
ProductBra
```

## Utility Functions

### Fock States

Convenience constructors for number states:

```@docs
FockKet
FockBra
```

## Examples

```julia
using QSymbolic

H, Hb = HilbertSpace(:H, 2)
Zb = Basis(H, :z)

# Create basis kets
up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Superposition
plus = (up + down) / √2

# Adjoint gives bra
up'  # → ⟨↑|

# Inner products
up' * up    # → 1
up' * down  # → 0

# Tensor products (order-independent)
H_A, Ba = HilbertSpace(:A, 2)
H_B, Bb = HilbertSpace(:B, 2)
ψ = Ket(Ba, :ψ)
ϕ = Ket(Bb, :ϕ)
ψ ⊗ ϕ == ϕ ⊗ ψ  # → true

# Symbolic indices
n = Sym(:n, :nonnegative, :integer)
F, Fb = FockSpace(:mode)
ket_n = Ket(Fb, n)  # |n⟩
```
