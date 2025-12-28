# Composite Systems

QSymbolic.jl supports tensor products for multi-particle quantum systems.

## Tensor Product of Spaces

Combine two Hilbert spaces:

```julia
using QSymbolic

H_A = HilbertSpace(:A, 2)
H_B = HilbertSpace(:B, 2)

H_AB = H_A ⊗ H_B  # CompositeSpace
```

## Product States

Create tensor products of kets:

```julia
ψ_A = BasisKet(H_A, :ψ)
ϕ_B = BasisKet(H_B, :ϕ)

product = ψ_A ⊗ ϕ_B  # |ψ⟩_A ⊗ |ϕ⟩_B
```

## Inner Products of Product States

Product state inner products factorize:

```math
(⟨ψ_A| ⊗ ⟨ϕ_B|)(|ψ'_A⟩ ⊗ |ϕ'_B⟩) = ⟨ψ_A|ψ'_A⟩ · ⟨ϕ_B|ϕ'_B⟩
```

```julia
product' * product  # → 1
```

## Composite Bases

When you have bases on each subsystem, they combine:

```julia
Za = Basis(H_A, :z)
Zb = Basis(H_B, :z)

ZaZb = Za ⊗ Zb  # CompositeBasis
```

## Factorized Transforms

The key feature: **transforms factorize automatically**.

If you have:
- Transform from `Xa → Za` (for subsystem A)
- Transform from `Xb → Zb` (for subsystem B)

Then `Xa⊗Xb → Za⊗Zb` is derived automatically:

```julia
# Define subsystem transforms
define_transform!(Xa, Za) do idx
    idx == :↑ ? (up_a + down_a)/√2 : (up_a - down_a)/√2
end
define_transform!(Xb, Zb) do idx
    idx == :↑ ? (up_b + down_b)/√2 : (up_b - down_b)/√2
end

# Product state in x⊗x basis
up_x_a = BasisKet(Xa, :↑)
up_x_b = BasisKet(Xb, :↑)
state_xx = up_x_a ⊗ up_x_b

# Transform to z⊗z basis - automatic!
state_zz = transform(state_xx, typeof(Za ⊗ Zb))
```

## Entangled States

Entangled states are superpositions of product states:

```julia
H = HilbertSpace(:qubit, 2)
Z = Basis(H, :z)

up = BasisKet(Z, :0)
down = BasisKet(Z, :1)

# Bell state |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
zero_zero = up ⊗ up
one_one = down ⊗ down
bell_plus = (zero_zero + one_one) / √2  # Not yet implemented for mixed products
```

!!! note
    Full support for entangled state arithmetic is in development.

## Explicit Composite Transforms

For non-factorizable transforms (e.g., to an entangled basis like Bell states), register explicitly:

```julia
BellBasis = Basis(H_A ⊗ H_B, :bell)  # hypothetical

define_transform!(Za ⊗ Zb, BellBasis) do ket
    # Custom entangled transform
end
```
