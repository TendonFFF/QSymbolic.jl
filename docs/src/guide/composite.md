# Composite Systems

QSymbolic.jl supports tensor products for multi-particle quantum systems, with automatic factorization of basis transforms.

## Tensor Product of Spaces

Combine two Hilbert spaces with the tensor product operator `⊗`:

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
ψ_A = BasisKet(H_A, :ψ)
ϕ_B = BasisKet(H_B, :ϕ)

product = ψ_A ⊗ ϕ_B  # |ψ⟩_A ⊗ |ϕ⟩_B (ProductKet)
```

This represents a state where subsystem A is in state |ψ⟩ and subsystem B is in state |ϕ⟩, with no correlations between them.

## Inner Products of Product States

Inner products of product states **factorize** according to:

```math
(⟨ψ_A| ⊗ ⟨ϕ_B|)(|ψ'_A⟩ ⊗ |ϕ'_B⟩) = ⟨ψ_A|ψ'_A⟩ \cdot ⟨ϕ_B|ϕ'_B⟩
```

```julia
product' * product  # → 1 × 1 = 1

# Different states
χ_A = BasisKet(H_A, :χ)
other = χ_A ⊗ ϕ_B

product' * other    # → ⟨ψ|χ⟩ × ⟨ϕ|ϕ⟩ = 0 × 1 = 0
```

## Composite Bases

When you create bases on each subsystem, they can be combined into a composite basis:

```julia
Za = Basis(H_A, :z)  # z-basis for qubit A
Zb = Basis(H_B, :z)  # z-basis for qubit B

ZaZb = Za ⊗ Zb       # CompositeBasis for the joint system
```

## Factorized Transforms

The key feature of QSymbolic.jl's composite system support: **transforms factorize automatically**.

If you define:
- Transform from `Xa → Za` (for subsystem A)
- Transform from `Xb → Zb` (for subsystem B)

Then the composite transform `Xa⊗Xb → Za⊗Zb` is derived automatically, without explicit registration!

### Example: Two-Qubit Basis Change

```julia
# Setup two qubits
H_A = HilbertSpace(:A, 2)
H_B = HilbertSpace(:B, 2)

# Define bases for each
Za, Xa = Basis(H_A, :z), Basis(H_A, :x)
Zb, Xb = Basis(H_B, :z), Basis(H_B, :x)

# Create states
up_a, down_a = BasisKet(Za, :↑), BasisKet(Za, :↓)
up_b, down_b = BasisKet(Zb, :↑), BasisKet(Zb, :↓)

# Define subsystem transforms
define_transform!(Xa, Za) do idx
    idx == :↑ ? (up_a + down_a)/√2 : (up_a - down_a)/√2
end
define_transform!(Xb, Zb) do idx
    idx == :↑ ? (up_b + down_b)/√2 : (up_b - down_b)/√2
end

# Now create a product state in x⊗x basis
up_x_a = BasisKet(Xa, :↑)
up_x_b = BasisKet(Xb, :↑)
state_xx = up_x_a ⊗ up_x_b  # |↑_x⟩_A ⊗ |↑_x⟩_B

# Transform to z⊗z basis - automatically factorized!
target_basis = typeof(Za ⊗ Zb)
state_zz = transform(state_xx, target_basis)
# Result: (|↑↑⟩ + |↑↓⟩ + |↓↑⟩ + |↓↓⟩)/2
```

The factorization uses the tensor product property:
```math
(U_A \otimes U_B)|ψ_A⟩|ψ_B⟩ = (U_A|ψ_A⟩) \otimes (U_B|ψ_B⟩)
```

## Entangled States

Entangled states are **superpositions** of product states that cannot be factored:

```julia
H = HilbertSpace(:qubit, 2)
Z = Basis(H, :z)

up = BasisKet(Z, :0)
down = BasisKet(Z, :1)

# Bell state |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
zero_zero = up ⊗ up      # |0⟩|0⟩
one_one = down ⊗ down    # |1⟩|1⟩

# Bell state as SumProductKet
bell_plus = (zero_zero + one_one) / √2
```

!!! warning "Entangled State Arithmetic"
    Arithmetic with `SumProductKet` (sums of product states) is partially supported. Inner products between `SumProductKet` states work correctly.

## Explicit Composite Transforms

For non-factorizable transforms (e.g., transforming to an entangled basis like Bell states), you can register transforms explicitly:

```julia
# Create a composite basis explicitly
CompBasis1 = Za ⊗ Zb

# For a hypothetical Bell basis, you would define:
# define_transform!(CompBasis1, BellBasis) do ket
#     # Custom transform to entangled basis
# end
```

## Checking Composite Transforms

```julia
# Check if factorized transform is available
has_transform(typeof(Xa ⊗ Xb), typeof(Za ⊗ Zb))  # → true (if both subsystem transforms exist)

# This works even without explicitly registering the composite transform!
```

## Tensor Product Operators

Operators on composite systems can be constructed as tensor products:

```julia
H_A = HilbertSpace(:A, 2)
H_B = HilbertSpace(:B, 2)
Za = Basis(H_A, :z)
Zb = Basis(H_B, :z)

up_a, down_a = BasisKet(Za, :↑), BasisKet(Za, :↓)
up_b, down_b = BasisKet(Zb, :↑), BasisKet(Zb, :↓)

# Single-qubit Pauli Z operators
σz_a = up_a * up_a' - down_a * down_a'
σz_b = up_b * up_b' - down_b * down_b'

# Tensor product: σz_A ⊗ σz_B
σz_ab = σz_a ⊗ σz_b

# Apply to product state
ψ = up_a ⊗ up_b
σz_ab * ψ  # → |↑↑⟩ (eigenvalue +1)
```

### Lifting Single-System Operators

Use `lift` to extend a single-system operator to a composite space with identity:

```julia
# Lift σz_a to act on joint system A⊗B
σz_a_full = lift(σz_a, Zb)  # equivalent to σz_a ⊗ 𝕀_B

# Using IdentityOp directly
σz_a_full = σz_a ⊗ IdentityOp(Zb)
σz_b_full = IdentityOp(Za) ⊗ σz_b
```

### Reordering Tensor Products

Use `reorder` to reorder tensor products to match a target basis ordering:

```julia
# If we have operator in order A⊗B but need B⊗A:
T_ab = σz_a ⊗ σz_b     # Acts on Za ⊗ Zb
T_ba = reorder(T_ab, (Zb, Za))  # Reordered for Zb ⊗ Za
```

### Swapping Subsystems

The `swap` function swaps two adjacent positions in a tensor product:

```julia
T = op1 ⊗ op2 ⊗ op3  # 3-system operator
T_swapped = swap(T, 1)  # Swap positions 1 and 2: becomes op2 ⊗ op1 ⊗ op3
```
