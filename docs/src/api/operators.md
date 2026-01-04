# Operators

Quantum operators transform states and are essential for describing observables and dynamics. QSymbolic.jl provides outer product operators, operator algebra, and function-defined operators.

## Overview

| Type | Description |
|:-----|:------------|
| `Operator` | Outer product `\|ψ⟩⟨ϕ\|` with coefficient |
| `SumOperator` | Sum of operators: Â + B̂ |
| `ScaledOperator` | Scalar times operator: α·Â |
| `OperatorProduct` | Product of operators: ÂB̂ |
| `FunctionOperator` | Operator defined by a function |
| `IdentityOp` | Identity operator 𝕀 |
| `TensorOperator` | Tensor product of operators: Â ⊗ B̂ |

## Abstract Type

```@docs
AbstractOperator
```

## Outer Product Operator

The primary way to build operators from states:

```@docs
Operator
```

## Operator Algebra Types

### Sum of Operators

```@docs
SumOperator
```

### Scaled Operator

```@docs
ScaledOperator
```

### Operator Product

```@docs
OperatorProduct
```

## Function-Based Operator

For operators with procedural definitions:

```@docs
FunctionOperator
AdjointFunctionOperator
```

## Identity Operator

```@docs
IdentityOp
```

## Tensor Product Operator

```@docs
TensorOperator
```

### Tensor Product Utilities

```@docs
lift
swap
reorder
partial_trace
```

## Symbolic Types

When operator application cannot be simplified:

```@docs
OpKet
OpBra
```

## Accessor Functions

```@docs
basis(::Operator)
```

## Examples

### Projectors and Ladder Operators

```julia
using QSymbolic

H = HilbertSpace(:spin, 2)
Hb = Basis(H, :default)
Zb = Basis(H, :z)
up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Projector
P_up = up * up'         # |↑⟩⟨↑|

# Ladder operators
σ_plus = up * down'     # |↑⟩⟨↓|
σ_minus = down * up'    # |↓⟩⟨↑|

# Apply
P_up * up       # → |↑⟩
σ_plus * down   # → |↑⟩
```

### Pauli Matrices

```julia
# Build from outer products
σx = up * down' + down * up'
σy = -im * (up * down') + im * (down * up')
σz = up * up' - down * down'

# Eigenvalue equations
σz * up    # → |↑⟩
σz * down  # → -|↓⟩
```

### Function Operator (Fock Space)

```julia
F = FockSpace(:mode)
Fb = Basis(F, :n)

# Annihilation operator
â = FunctionOperator(:â, Fb) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end

# Creation operator  
â† = FunctionOperator(:â†, Fb) do ket
    n = parse(Int, string(ket.index))
    √(n + 1) * Ket(Fb, n + 1)
end
```

### Tensor Product Operators

```julia
# Two-qubit system
H1 = HilbertSpace(:qubit1, 2)
H1b = Basis(H1, :default)
H2 = HilbertSpace(:qubit2, 2)
H2b = Basis(H2, :default)
B1 = Basis(H1, :z)
B2 = Basis(H2, :z)

up1 = Ket(B1, :↑)
down1 = Ket(B1, :↓)
up2 = Ket(B2, :↑)
down2 = Ket(B2, :↓)

# Single-qubit operators
σz1 = up1 * up1' - down1 * down1'
σz2 = up2 * up2' - down2 * down2'

# Tensor product
σz1_σz2 = σz1 ⊗ σz2  # σz ⊗ σz

# Lift operator to composite space with identity
σz1_full = σz1 ⊗ IdentityOp(B2)  # σz ⊗ 𝕀

# Using lift function
σz1_lifted = lift(σz1, B2)  # equivalent to σz1 ⊗ 𝕀(B2)

# Reorder tensor product to match target basis order
T12 = σz1 ⊗ σz2
T21 = reorder(T12, (B2, B1))  # reorder to B2⊗B1
```
