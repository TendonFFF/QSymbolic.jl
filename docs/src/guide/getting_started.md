# Getting Started

This guide introduces the basic concepts of QSymbolic.jl.

## Hilbert Spaces

A Hilbert space is the fundamental arena for quantum mechanics. Create one with:

```julia
using QSymbolic

# Finite-dimensional (e.g., qubit)
H = HilbertSpace(:qubit, 2)

# Infinite-dimensional (e.g., harmonic oscillator)
F = FockSpace(:oscillator)
```

## States: Kets and Bras

Quantum states are represented as kets ``|ψ⟩`` and their duals as bras ``⟨ψ|``:

```julia
H = HilbertSpace(:H, 2)

# Create kets
ψ = BasisKet(H, :ψ)
ϕ = BasisKet(H, :ϕ)

# Bras via adjoint
ψ'  # ⟨ψ|

# Inner products
ψ' * ψ  # → 1 (same state)
ψ' * ϕ  # → 0 (orthogonal by default)
```

## Linear Combinations

Build superpositions with arithmetic:

```julia
# Scalar multiplication
2 * ψ       # 2|ψ⟩

# Addition/subtraction
ψ + ϕ       # |ψ⟩ + |ϕ⟩
ψ - ϕ       # |ψ⟩ - |ϕ⟩

# Normalization
(ψ + ϕ) / √2
```

## Inner Products

Inner products follow Dirac notation:

```julia
# Single states
ψ' * ψ  # → 1

# Superpositions
(ψ + ϕ)' * (ψ + ϕ)  # → 2 (since ⟨ψ|ψ⟩ + ⟨ϕ|ϕ⟩ = 1 + 1)
```

## Next Steps

- Learn about [Basis Transforms](@ref) for changing representations
- Explore [Composite Systems](@ref) for multi-particle states
