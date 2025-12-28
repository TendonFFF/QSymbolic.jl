# Bases

Bases provide the reference frame for expressing quantum states. QSymbolic.jl supports explicit named bases with automatic orthonormality.

## Overview

| Type | Description |
|:-----|:------------|
| `Basis` | A named orthonormal basis for a Hilbert space |
| `DefaultBasis` | Implicit basis used when none is specified |
| `CompositeBasis` | Tensor product of two bases |

## Types

```@docs
AbstractBasis
Basis
DefaultBasis
CompositeBasis
```

## Accessor Functions

```@docs
space
basisname
basis1
basis2
```

## Tensor Product

Use `⊗` to combine bases:

```@docs
⊗(::AbstractBasis, ::AbstractBasis)
```

## Examples

```julia
using QSymbolic

H = HilbertSpace(:spin, 2)

# Create named bases
Zb = Basis(H, :z)
Xb = Basis(H, :x)

# Get the space a basis is defined on
space(Zb)  # HilbertSpace(:spin, 2)

# Get the basis name
basisname(Zb)  # :z

# Composite basis
H_A, H_B = HilbertSpace(:A, 2), HilbertSpace(:B, 2)
ZaZb = Basis(H_A, :z) ⊗ Basis(H_B, :z)
```
