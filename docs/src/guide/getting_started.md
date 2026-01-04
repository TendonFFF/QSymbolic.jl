# Getting Started

This guide introduces the basic concepts of QSymbolic.jl for symbolic quantum mechanics.

## Hilbert Spaces and Bases

A Hilbert space is the fundamental mathematical arena for quantum mechanics. In QSymbolic.jl, you create spaces with the `HilbertSpace` constructor:

```julia
using QSymbolic

# Finite-dimensional space (e.g., qubit with 2 dimensions)
H = HilbertSpace(:qubit, 2)

# Infinite-dimensional space (e.g., harmonic oscillator)
F = HilbertSpace(:oscillator)

# Convenience alias for Fock spaces
F = FockSpace(:oscillator)  # equivalent to HilbertSpace(:oscillator)
```

### Creating Bases

Every ket lives in a **basis**. Create named bases with `Basis(space, name)`:

```julia
H = HilbertSpace(:spin, 2)

# Create named bases
Zb = Basis(H, :z)   # spin-z eigenbasis
Xb = Basis(H, :x)   # spin-x eigenbasis
```

### Convenient Destructuring

`HilbertSpace` supports destructuring to get both a space and a default basis in one line:

```julia
# Get space AND default basis
H, Hb = HilbertSpace(:H, 2)

# Equivalent to:
# H = HilbertSpace(:H, 2)
# Hb = Basis(H, :default)

# Same works for FockSpace
F, Fb = FockSpace(:mode)
```

This pattern is convenient and used throughout the documentation.

## States: Kets and Bras

Quantum states are represented in Dirac notation. A **ket** `|ψ⟩` represents a quantum state, and its dual **bra** `⟨ψ|` is the conjugate transpose.

### Creating Kets

```julia
H, Hb = HilbertSpace(:H, 2)

# Create kets with symbolic labels
ψ = Ket(Hb, :ψ)
ϕ = Ket(Hb, :ϕ)

# Display uses Dirac notation
ψ   # → |ψ⟩

# Get the bra (conjugate transpose) via adjoint
ψ'  # → ⟨ψ|
```

### Numeric Indices (Fock States)

For number states, use integer indices:

```julia
F, Fb = FockSpace(:mode)

# Number states |n⟩
n0 = Ket(Fb, 0)   # Ground state |0⟩
n1 = Ket(Fb, 1)   # First excited |1⟩
n2 = Ket(Fb, 2)   # Second excited |2⟩

# Orthonormality
n0' * n0  # → 1
n0' * n1  # → 0
```

### Symbolic Indices

Kets can have symbolic indices for general representations:

```julia
F, Fb = FockSpace(:mode)

# Create symbolic index
n = Sym(:n, :nonnegative, :integer)

# Ket with symbolic index
ket_n = Ket(Fb, n)   # |n⟩

# Inner products give KroneckerDelta
m = Sym(:m, :nonnegative, :integer)
Ket(Fb, m)' * ket_n  # → δ(m,n)
```

!!! note "Orthonormality"
    States in the same basis with different labels are automatically orthogonal. This reflects the standard assumption that basis states form an orthonormal set.

## Inner Products

Inner products follow quantum mechanics conventions with `⟨bra|ket⟩` notation:

```julia
H, Hb = HilbertSpace(:H, 2)
ψ = Ket(Hb, :ψ)
ϕ = Ket(Hb, :ϕ)

# Orthonormal basis states
ψ' * ψ  # → 1 (normalized)
ψ' * ϕ  # → 0 (orthogonal - different labels)

# Inner product of superpositions expands linearly
# ⟨ψ + ϕ|ψ + ϕ⟩ = ⟨ψ|ψ⟩ + ⟨ψ|ϕ⟩ + ⟨ϕ|ψ⟩ + ⟨ϕ|ϕ⟩ = 1 + 0 + 0 + 1 = 2
(ψ + ϕ)' * (ψ + ϕ)  # → 2
```

## Linear Combinations (Superpositions)

Quantum superpositions are built using standard arithmetic:

```julia
H, Hb = HilbertSpace(:H, 2)
ψ = Ket(Hb, :ψ)
ϕ = Ket(Hb, :ϕ)

# Scalar multiplication
2 * ψ           # → 2|ψ⟩

# Addition creates superpositions
ψ + ϕ           # → |ψ⟩ + |ϕ⟩

# Subtraction
ψ - ϕ           # → |ψ⟩ - |ϕ⟩

# Division by scalars (for normalization)
(ψ + ϕ) / √2    # → (|ψ⟩ + |ϕ⟩)/√2

# Rational division also works
ψ // 2          # → (1//2)|ψ⟩
```

### Type Promotion

Operations automatically create appropriate types:

| Operation | Input | Result Type |
|:----------|:------|:------------|
| `α * ψ` | Scalar × Ket | `WeightedKet` |
| `ψ + ϕ` | Ket + Ket | `SumKet` |
| `ψ ⊗ ϕ` | Ket ⊗ Ket | `ProductKet` |

## Symbolic Scalars

QSymbolic.jl uses [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) for symbolic computation. Create symbolic variables with `Sym`:

```julia
# Basic symbolic variable
n = Sym(:n)

# With type assumptions
θ = Sym(:θ, :real)
k = Sym(:k, :positive, :integer)

# Using Symbolics.jl macros (re-exported)
@variables x y z
@syms α β
```

### Assumptions

Assumptions affect simplification and conjugation:

```julia
# Complex (default) - conjugate applies
z = Sym(:z)
z'  # → conj(z)

# Real - self-conjugate
r = Sym(:r, :real)
r'  # → r

# Available: :real, :positive, :negative, :nonnegative, :integer
```

### Symbolic State Coefficients

```julia
H, Hb = HilbertSpace(:H, 2)
ψ = Ket(Hb, :ψ)
ϕ = Ket(Hb, :ϕ)

# Symbolic amplitudes
α = Sym(:α)
β = Sym(:β)

# Superposition with symbolic coefficients
state = α * ψ + β * ϕ  # → α|ψ⟩ + β|ϕ⟩
```

## Outer Products and Operators

Create operators via outer products `|ψ⟩⟨ϕ|`:

```julia
H, Hb = HilbertSpace(:spin, 2)
Zb = Basis(H, :z)

up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Projector onto |↑⟩
P_up = up * up'      # |↑⟩⟨↑|

# Apply to states
P_up * up    # → |↑⟩
P_up * down  # → 0

# Ladder operator
σ_plus = up * down'  # |↑⟩⟨↓|
σ_plus * down  # → |↑⟩

# Pauli Z
σz = up * up' - down * down'
σz * up    # → |↑⟩
σz * down  # → -|↓⟩
```

## Quick Example: Spin-1/2 System

```julia
using QSymbolic

# Create 2D Hilbert space with z-basis
H = HilbertSpace(:spin, 2)
Zb = Basis(H, :z)

# Basis states
up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Superposition (equal superposition)
ψ = (up + down) / √2

# Check normalization
ψ' * ψ  # → 1

# Pauli operators
σz = up * up' - down * down'
σx = up * down' + down * up'

# Apply operators
σz * up     # → |↑⟩
σz * down   # → -|↓⟩
σx * up     # → |↓⟩
σx * down   # → |↑⟩

# Measurement probability
# P(↑) = |⟨↑|ψ⟩|² 
up' * ψ  # → 1/√2, so P(↑) = 1/2
```

## Next Steps

- [Basis Transforms](transforms.md) - Changing between representations (z ↔ x basis)
- [Composite Systems](@ref) - Multi-particle states and tensor products
- [Operators](operators.md) - Function operators and operator algebra
- [Symbolic Scalars](symbolic.md) - Advanced symbolic computation
- [Custom Contraction Rules](@ref) - Non-orthonormal bases
