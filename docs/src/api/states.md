# States

Quantum states are represented as kets `|ψ⟩` and their duals as bras `⟨ψ|`. QSymbolic.jl provides types for basis states, superpositions, and product states.

## Overview

| Type | Description |
|:-----|:------------|
| `BasisKet` / `BasisBra` | Single basis state |
| `weightedKet` / `weightedBra` | Scalar × basis state |
| `sumKet` / `sumBra` | Linear combination of basis states |
| `ProductKet` / `ProductBra` | Tensor product of states |
| `SumProductKet` / `SumProductBra` | Superposition of product states |
| `InnerProduct` | Symbolic (unevaluated) inner product |

## Single-System States

### Abstract Types

```@docs
AbstractKet
AbstractBra
```

### Basis States

```@docs
BasisKet
BasisBra
```

### Weighted States

```@docs
weightedKet
weightedBra
```

### Superpositions

```@docs
sumKet
sumBra
```

## Composite States

For multi-particle systems:

```@docs
ProductKet
ProductBra
SumProductKet
SumProductBra
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
up = BasisKet(Zb, :↑)
down = BasisKet(Zb, :↓)

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
