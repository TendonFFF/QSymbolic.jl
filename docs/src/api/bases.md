# Bases

Bases provide the reference frame for expressing quantum states. QSymbolic.jl supports explicit named bases with automatic orthonormality.

## Overview

| Type | Description |
|:-----|:------------|
| `Basis` | A named orthonormal basis for a Hilbert space |
| `DefaultBasis` | Implicit basis from `HilbertSpace` destructuring |
| `CompositeBasis` | Tensor product of bases |

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
bases
```

## Tensor Product

Use `⊗` to combine bases:

```@docs
⊗(::AbstractBasis, ::AbstractBasis)
```

## Examples

```julia
using QSymbolic

# Create space with destructuring (gets default basis)
H, Hb = HilbertSpace(:spin, 2)

# Create named bases
Zb = Basis(H, :z)
Xb = Basis(H, :x)

# Get the space a basis is defined on
space(Zb)  # → HilbertSpace{(:spin,), (2,)}

# Get the basis name
basisname(Zb)  # → :z

# Composite basis
H_A = HilbertSpace(:A, 2)
H_B = HilbertSpace(:B, 2)
Za = Basis(H_A, :z)
Zb = Basis(H_B, :z)
ZaZb = Za ⊗ Zb  # CompositeBasis

# Access components
basis1(ZaZb)  # → typeof(Za)
basis2(ZaZb)  # → typeof(Zb)

# Three-way composite
H_C = HilbertSpace(:C, 2)
Zc = Basis(H_C, :z)
ZaZbZc = Za ⊗ Zb ⊗ Zc
bases(ZaZbZc)  # → (Za, Zb, Zc) types
```
