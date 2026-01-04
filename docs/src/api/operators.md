# Operators

Quantum operators transform states and are essential for describing observables and dynamics. QSymbolic.jl provides outer product operators, operator algebra, and function-defined operators.

## Overview

| Type | Description |
|:-----|:------------|
| `Outer` | Single outer product `\|ψ⟩⟨ϕ\|` |
| `Operator` | Sum of weighted outer products |
| `Identity` | Identity operator 𝕀 |
| `FunctionOperator` | Operator defined by a function |

## Abstract Type

```@docs
AbstractOperator
```

## Outer Product Operator

The primary way to build operators from states:

```@docs
Outer
```

## Operator Container

Sum of weighted outer products:

```@docs
Operator
```

## Identity Operator

```@docs
Identity
```

## Function-Based Operator

For operators with procedural definitions (e.g., Fock space ladder operators):

```@docs
FunctionOperator
```

### FunctionOperator Syntax

```julia
# Basic usage
action(ket::Ket{B}) where B = ...  # returns AbstractKet or Number
op = FunctionOperator(action, basis, name=:op)

# With adjoint action
op = FunctionOperator(action, basis, adjoint_action=adj_action, name=:op)

# Do-block syntax
op = FunctionOperator(basis, name=:op) do ket
    # action on ket
end
```

## Accessor Functions

```@docs
space(::AbstractOperator)
```

## Examples

### Projectors and Ladder Operators

```julia
using QSymbolic

H, Hb = HilbertSpace(:spin, 2)
Zb = Basis(H, :z)
up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Projector
P_up = up * up'  # |↑⟩⟨↑|
P_up * up   # → |↑⟩
P_up * down # → 0

# Ladder operator
σ_plus = up * down'  # |↑⟩⟨↓|
σ_plus * down  # → |↑⟩

# Operator sum
σz = up * up' - down * down'
```

### Fock Space Operators

```julia
using QSymbolic

F, Fb = FockSpace(:mode)

# Annihilation operator: â|n⟩ = √n |n-1⟩
annihilate(ket::Ket{B}) where B = √(ket.index) * Ket{B}(ket.index - 1)
create(ket::Ket{B}) where B = √(ket.index + 1) * Ket{B}(ket.index + 1)

â = FunctionOperator(annihilate, Fb, adjoint_action=create, name=:â)

# Apply to symbolic Fock state
n = Sym(:n, :nonnegative, :integer)
ket_n = Ket(Fb, n)
â * ket_n   # → √n |n-1⟩
â' * ket_n  # → √(n+1) |n+1⟩
```
