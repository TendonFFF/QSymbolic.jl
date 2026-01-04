# States

Quantum states are represented as kets `|ψ⟩` and their duals as bras `⟨ψ|`. QSymbolic.jl provides types for basis states, superpositions, and product states.

## Overview

| Type | Description |
|:-----|:------------|
| `Ket` / `Bra` | Single basis state |
| `WeightedKet` / `WeightedBra` | Scalar × basis state |
| `SumKet` / `SumBra` | Linear combination of basis states |
| `ProductKet` / `ProductBra` | Tensor product of states |
| `SumKet` / `SumBra` | Superposition of product states |
| `InnerProduct` | Symbolic (unevaluated) inner product |

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

For multi-particle systems:

```@docs
ProductKet
ProductBra
SumKet
SumBra
```

## Symbolic Types

When inner products cannot be evaluated (e.g., cross-basis without a transform):

```@docs
InnerProduct
```

## Utility Functions

### Accessing Basis and Space

```@docs
basis
```

### Fock States

Convenience constructors for number states in infinite-dimensional spaces:

```@docs
FockKet
FockBra
```

### Validation

```@docs
check_basis
check_space
```

## Examples

```julia
using QSymbolic

H = HilbertSpace(:H, 2)
Zb = Basis(H, :z)

# Create basis kets
up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Superposition
plus = (up + down) / √2

# Get the basis
basis(up)  # Basis{HilbertSpace{(:H,), (2,)}, :z}

# Adjoint gives bra
up'  # ⟨↑|

# Inner products
up' * up    # 1
up' * down  # 0
```
