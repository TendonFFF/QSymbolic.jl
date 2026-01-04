# Composite Systems

QSymbolic.jl supports tensor products for multi-particle quantum systems, with automatic factorization of basis transforms.

## Tensor Product of Spaces

Combine Hilbert spaces with the tensor product operator `⊗`:

```julia
using QSymbolic

H_A = HilbertSpace(:A, 2)  # First qubit
H_B = HilbertSpace(:B, 2)  # Second qubit

H_AB = H_A ⊗ H_B  # 4-dimensional composite space
```

The resulting `CompositeSpace` represents the joint system.

## Product States

Create tensor products of kets for separable (non-entangled) states:

```julia
H_A, Ba = HilbertSpace(:A, 2)
H_B, Bb = HilbertSpace(:B, 2)

ψ_A = Ket(Ba, :ψ)
ϕ_B = Ket(Bb, :ϕ)

product = ψ_A ⊗ ϕ_B  # |ψ⟩_A ⊗ |ϕ⟩_B (ProductKet)
```

This represents a state where subsystem A is in state |ψ⟩ and subsystem B is in state |ϕ⟩, with no correlations between them.

## Order-Independent (Bosonic) Behavior

!!! note "ProductKet Ordering"
    `ProductKet`s are **order-independent** (bosonic/symmetric). The kets are canonically ordered by basis, so `ψ_A ⊗ ϕ_B == ϕ_B ⊗ ψ_A` when they have different bases. This is the default behavior.

```julia
H_A, Ba = HilbertSpace(:A, 2)
H_B, Bb = HilbertSpace(:B, 2)

ψ = Ket(Ba, :ψ)
ϕ = Ket(Bb, :ϕ)

# Order doesn't matter for different bases
ψ ⊗ ϕ == ϕ ⊗ ψ  # → true

# Three-way tensor products also canonicalize
H_C, Bc = HilbertSpace(:C, 2)
χ = Ket(Bc, :χ)

ψ ⊗ ϕ ⊗ χ == χ ⊗ ψ ⊗ ϕ  # → true
```

!!! warning "Future Feature"
    Fermionic (anti-symmetric) tensor products will be added in a future update.

## Inner Products of Product States

Inner products of product states **factorize**:

```math
(⟨ψ_A| ⊗ ⟨ϕ_B|)(|ψ'_A⟩ ⊗ |ϕ'_B⟩) = ⟨ψ_A|ψ'_A⟩ \cdot ⟨ϕ_B|ϕ'_B⟩
```

```julia
product' * product  # → 1 × 1 = 1

# Different states
χ_A = Ket(Ba, :χ)
other = χ_A ⊗ ϕ_B

product' * other    # → ⟨ψ|χ⟩ × ⟨ϕ|ϕ⟩ = 0 × 1 = 0
```

## Multi-Index Kets

For composite systems, you can also use **multi-index kets** with a composite basis:

```julia
S_cavity, B_cavity = FockSpace(:cavity)
S_dot, B_dot = HilbertSpace(:dot, 3)

# Create composite space and basis
S_system = S_cavity ⊗ S_dot
B_composite = Basis(S_system, :dressed)

# Multi-index ket: |n, σ⟩
ket = Ket(B_composite, (Sym(:n), Sym(:σ)))

# Inner products check all indices
ket' * ket  # → 1
```

This is useful for dressed states like Jaynes-Cummings eigenstates. See [Custom Contraction Rules](@ref) for defining custom inner product behavior.

## Composite Bases

Combine bases into a composite basis:

```julia
H_A = HilbertSpace(:A, 2)
H_B = HilbertSpace(:B, 2)

Za = Basis(H_A, :z)
Zb = Basis(H_B, :z)

ZaZb = Za ⊗ Zb  # CompositeBasis for the joint system
```

## Factorized Transforms

**Transforms factorize automatically**. If you define:
- Transform from `Xa → Za` (for subsystem A)
- Transform from `Xb → Zb` (for subsystem B)

Then the composite transform `Xa⊗Xb → Za⊗Zb` is derived automatically!

```julia
H_A = HilbertSpace(:A, 2)
H_B = HilbertSpace(:B, 2)

# Define bases
Za, Xa = Basis(H_A, :z), Basis(H_A, :x)
Zb, Xb = Basis(H_B, :z), Basis(H_B, :x)

# Create states
up_a, down_a = Ket(Za, :↑), Ket(Za, :↓)
up_b, down_b = Ket(Zb, :↑), Ket(Zb, :↓)

# Define subsystem transforms
define_transform!(Xa, Za) do idx
    idx == :↑ ? (up_a + down_a)/√2 : (up_a - down_a)/√2
end
define_transform!(Xb, Zb) do idx
    idx == :↑ ? (up_b + down_b)/√2 : (up_b - down_b)/√2
end

# Product state in x⊗x basis
up_x_a = Ket(Xa, :↑)
up_x_b = Ket(Xb, :↑)
state_xx = up_x_a ⊗ up_x_b

# Transform to z⊗z - automatically factorized!
target_basis = typeof(Za ⊗ Zb)
state_zz = transform(state_xx, target_basis)
# Result: (|↑↑⟩ + |↑↓⟩ + |↓↑⟩ + |↓↓⟩)/2
```

## Entangled States

Entangled states are **superpositions** of product states:

```julia
H = HilbertSpace(:qubit, 2)
Z = Basis(H, :z)

up = Ket(Z, :0)
down = Ket(Z, :1)

# Bell state |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
zero_zero = up ⊗ up
one_one = down ⊗ down

bell_plus = (zero_zero + one_one) / √2
```

## Operators on Composite Systems

### Tensor Product of Operators

```julia
H_A, Ba = HilbertSpace(:A, 2)
H_B, Bb = HilbertSpace(:B, 2)
Za = Basis(H_A, :z)
Zb = Basis(H_B, :z)

up_a, down_a = Ket(Za, :↑), Ket(Za, :↓)
up_b, down_b = Ket(Zb, :↑), Ket(Zb, :↓)

# Single-qubit operators
σz_a = up_a * up_a' - down_a * down_a'
σz_b = up_b * up_b' - down_b * down_b'

# Tensor product: σz_A ⊗ σz_B
σz_ab = σz_a ⊗ σz_b

# Apply to product state
ψ = up_a ⊗ up_b
σz_ab * ψ  # → |↑↑⟩ (eigenvalue +1)
```

### Single-System Operators on Product States

Operators automatically act on their matching subsystem:

```julia
S_cavity, B_cavity = FockSpace(:cavity)
S_dot, B_dot = HilbertSpace(:dot, 3)

# Annihilation operator on cavity
annihilate(ket::Ket{B}) where B = √(ket.index) * Ket{B}(ket.index - 1)
a = FunctionOperator(annihilate, B_cavity, name=:a)

# Apply to product state
d = Sym(:d)
product_state = Ket(B_cavity, d) ⊗ Ket(B_dot, :g)

a * product_state  # → √d |d-1⟩ ⊗ |g⟩
```

## See Also

- [Basis Transforms](transforms.md) - Cross-basis computations
- [Operators](operators.md) - Operator algebra
- [Custom Contraction Rules](@ref) - Multi-index inner products
