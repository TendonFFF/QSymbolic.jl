# Spaces

Hilbert spaces form the mathematical foundation of quantum mechanics. QSymbolic.jl provides types to represent both simple and composite spaces.

## Overview

| Type | Description |
|:-----|:------------|
| `HilbertSpace` | A single named Hilbert space (finite or infinite-dimensional) |
| `CompositeSpace` | Tensor product of two or more spaces |
| `FockSpace` | Convenience constructor for infinite-dimensional spaces |

## Types

```@docs
AbstractSpace
HilbertSpace
CompositeSpace
FockSpace
```

## Tensor Product

Use `⊗` (typed as `\otimes<tab>`) to create composite spaces:

```@docs
⊗(::AbstractSpace, ::AbstractSpace)
```

## Examples

```julia
using QSymbolic

# Single qubit
H = HilbertSpace(:qubit, 2)

# Harmonic oscillator (infinite-dimensional)  
F = FockSpace(:oscillator)

# Two-qubit system
H_AB = HilbertSpace(:A, 2) ⊗ HilbertSpace(:B, 2)
```
