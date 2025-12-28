# Basis Transforms

QSymbolic.jl supports explicit bases and transformations between them.

## Why Explicit Bases?

In quantum mechanics, the same state can be expressed in different bases. For example, a spin-1/2 particle's "spin-up" state in the z-basis can be written as a superposition in the x-basis:

```math
|↑_z⟩ = \frac{1}{\sqrt{2}}(|↑_x⟩ + |↓_x⟩)
```

QSymbolic.jl makes this explicit with named bases.

## Creating Bases

```julia
using QSymbolic

H = HilbertSpace(:spin, 2)

# Define two different bases for the same space
Zb = Basis(H, :z)  # spin-z eigenbasis
Xb = Basis(H, :x)  # spin-x eigenbasis

# Create states in specific bases
up_z = BasisKet(Zb, :↑)
down_z = BasisKet(Zb, :↓)
up_x = BasisKet(Xb, :↑)
```

## Orthonormality

States in the **same basis** are automatically orthonormal:

```julia
up_z' * up_z    # → 1
up_z' * down_z  # → 0
```

## Cross-Basis Inner Products

Without a transform, cross-basis inner products return a symbolic result:

```julia
up_z' * up_x  # → InnerProduct (symbolic)
```

## Registering Transforms

Define how one basis transforms to another:

```julia
define_transform!(Xb, Zb) do idx
    if idx == :↑
        (up_z + down_z) / √2
    else
        (up_z - down_z) / √2
    end
end
```

Now cross-basis inner products compute:

```julia
up_z' * up_x   # → 1/√2
down_z' * up_x # → 1/√2
```

## Transform API

```julia
# Check if transform exists
has_transform(typeof(Xb), typeof(Zb))  # → true

# Apply transform explicitly
transform(up_x, typeof(Zb))  # → (|↑⟩ + |↓⟩)/√2

# Clear all transforms
clear_transforms!()
```

## Inverse Transforms

Transforms are directional. If you need both directions, register both:

```julia
# X → Z
define_transform!(Xb, Zb) do idx
    idx == :↑ ? (up_z + down_z)/√2 : (up_z - down_z)/√2
end

# Z → X
define_transform!(Zb, Xb) do idx
    idx == :↑ ? (up_x + down_x)/√2 : (up_x - down_x)/√2
end
```
