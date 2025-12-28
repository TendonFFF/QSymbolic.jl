# Getting Started

This guide introduces the basic concepts of QSymbolic.jl for symbolic quantum mechanics.

## Hilbert Spaces

A Hilbert space is the fundamental mathematical arena for quantum mechanics. In QSymbolic.jl, you create spaces with the `HilbertSpace` constructor:

```julia
using QSymbolic

# Finite-dimensional space (e.g., qubit with 2 dimensions)
H = HilbertSpace(:qubit, 2)

# Infinite-dimensional space (e.g., harmonic oscillator)
F = HilbertSpace(:oscillator)

# Convenience alias for Fock spaces
F = FockSpace(:oscillator)  # equivalent to above
```

The first argument is a symbolic name, and the optional second argument specifies the dimension.

## States: Kets and Bras

Quantum states are represented in Dirac notation. A **ket** ``|ψ⟩`` represents a quantum state, and its dual **bra** ``⟨ψ|`` is the conjugate transpose:

```julia
H = HilbertSpace(:H, 2)

# Create kets with symbolic labels
ψ = BasisKet(H, :ψ)
ϕ = BasisKet(H, :ϕ)

# Display uses Dirac notation
ψ   # |ψ⟩

# Get the bra (conjugate transpose) via adjoint
ψ'  # ⟨ψ|

# Inner products follow ⟨bra|ket⟩ notation
ψ' * ψ  # → 1 (same state, normalized)
ψ' * ϕ  # → 0 (different labels are orthogonal by default)
```

!!! note "Orthonormality"
    States in the same basis with different labels are automatically orthogonal. This reflects the standard assumption that basis states form an orthonormal set.

## Linear Combinations (Superpositions)

Quantum superpositions are built using standard arithmetic:

```julia
# Scalar multiplication
2 * ψ           # 2|ψ⟩
(1/√2) * ψ      # (1/√2)|ψ⟩

# Addition creates superpositions
ψ + ϕ           # |ψ⟩ + |ϕ⟩

# Subtraction
ψ - ϕ           # |ψ⟩ - |ϕ⟩

# Division by scalars (for normalization)
(ψ + ϕ) / √2    # (|ψ⟩ + |ϕ⟩)/√2

# Rational division also works
ψ // 2          # (1//2)|ψ⟩
```

## Inner Products

Inner products follow quantum mechanics conventions:

```julia
H = HilbertSpace(:H, 2)
ψ = BasisKet(H, :ψ)
ϕ = BasisKet(H, :ϕ)

# Orthonormal basis states
ψ' * ψ  # → 1 (normalized)
ψ' * ϕ  # → 0 (orthogonal)

# Inner product of superpositions expands linearly
# ⟨ψ + ϕ|ψ + ϕ⟩ = ⟨ψ|ψ⟩ + ⟨ψ|ϕ⟩ + ⟨ϕ|ψ⟩ + ⟨ϕ|ϕ⟩ = 1 + 0 + 0 + 1 = 2
(ψ + ϕ)' * (ψ + ϕ)  # → 2
```

## Fock States

For quantum optics and harmonic oscillators, use Fock states with integer labels:

```julia
F = FockSpace(:mode)

# Number states |n⟩
n0 = FockKet(F, 0)  # Ground state |0⟩
n1 = FockKet(F, 1)  # First excited state |1⟩
n2 = FockKet(F, 2)  # Second excited state |2⟩

# Orthonormality
n0' * n0  # → 1
n0' * n1  # → 0
```

## Symbolic Scalars with Assumptions

QSymbolic.jl provides symbolic scalars that can carry type assumptions, which affects how operations like conjugation behave:

```julia
# Basic symbolic variable (no assumptions - treated as complex)
z = Sym(:z)
z'  # → z* (conjugate applied)

# Real variable - conjugate is identity
r = Sym(:r, :real)
r'  # → r (no conjugate for real numbers)

# Positive integer - implies real
n = Sym(:n, :positive, :integer)
is_real(n)       # → true
is_positive(n)   # → true
is_integer(n)    # → true
n'               # → n (adjoint of real is itself)

# Keyword syntax also supported
k = Sym(:k, positive=true, integer=true)
```

Available assumptions:
- `:real` - real number (conj(x) = x)
- `:positive` - strictly positive (implies :real and :nonnegative)
- `:negative` - strictly negative (implies :real)
- `:nonnegative` - non-negative (implies :real)
- `:integer` - integer value (implies :real)

This is particularly useful when working with quantum number indices:

```julia
# Symbolic quantum number
n = Sym(:n, :nonnegative, :integer)

# Expression involving the quantum number
expr = √(n + 1)  # √(n + 1) - properly real-valued
```

## Next Steps

- Learn about [Basis Transforms](@ref) for changing between different representations (e.g., position ↔ momentum)
- Explore [Composite Systems](@ref) for multi-particle states and tensor products
