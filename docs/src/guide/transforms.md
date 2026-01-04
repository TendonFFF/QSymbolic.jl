# Basis Transforms

QSymbolic.jl supports explicit bases and user-defined transformations between them. This is essential for quantum mechanics, where the same physical state can be expressed in different representations.

## Why Explicit Bases?

In quantum mechanics, choosing a basis is like choosing a coordinate system. The same state looks different in different bases. For example, a spin-1/2 particle's "spin-up" state in the z-basis is a superposition in the x-basis:

```math
|↑_z⟩ = \frac{1}{\sqrt{2}}(|↑_x⟩ + |↓_x⟩)
```

QSymbolic.jl makes bases explicit, so you can:
- Track which basis states live in
- Define transforms between bases
- Automatically compute cross-basis inner products

## Creating Bases

Use `Basis(space, name)` to create a named orthonormal basis:

```julia
using QSymbolic

H = HilbertSpace(:spin, 2)
Hb = Basis(H, :default)

# Define two different bases for the same space
Zb = Basis(H, :z)  # spin-z eigenbasis
Xb = Basis(H, :x)  # spin-x eigenbasis

# Create states in specific bases
up_z = Ket(Zb, :↑)
down_z = Ket(Zb, :↓)
up_x = Ket(Xb, :↑)
down_x = Ket(Xb, :↓)
```

!!! tip "Default Basis"
    When you use `Ket(space, label)` directly (without an explicit basis), QSymbolic uses an implicit `DefaultBasis`. This is convenient for simple calculations where you don't need multiple bases.

## Orthonormality Within a Basis

States in the **same basis** are automatically orthonormal:

```julia
up_z' * up_z    # → 1 (normalized)
up_z' * down_z  # → 0 (orthogonal)

up_x' * up_x    # → 1
up_x' * down_x  # → 0
```

## Cross-Basis Inner Products

What happens when you compute ``⟨↑_z|↑_x⟩``? Without knowing the relationship between bases, this is undefined. QSymbolic returns a symbolic result:

```julia
up_z' * up_x  # → ⟨↑|↑⟩ (InnerProduct - symbolic, unevaluated)
```

## Registering Transforms

To compute cross-basis inner products, tell QSymbolic how to transform between bases using `define_transform!`:

```julia
# Define transform from x-basis to z-basis
# The function takes an index and returns the state in the target basis
define_transform!(Xb, Zb) do idx
    if idx == :↑
        (up_z + down_z) / √2    # |↑_x⟩ = (|↑_z⟩ + |↓_z⟩)/√2
    else
        (up_z - down_z) / √2    # |↓_x⟩ = (|↑_z⟩ - |↓_z⟩)/√2
    end
end
```

Now cross-basis inner products are computed automatically:

```julia
up_z' * up_x    # → 0.7071... (1/√2)
down_z' * up_x  # → 0.7071... (1/√2)
up_z' * down_x  # → 0.7071... (1/√2)
down_z' * down_x # → -0.7071... (-1/√2)
```

## Transform API

### Check if a Transform Exists

```julia
has_transform(typeof(Xb), typeof(Zb))  # → true
has_transform(typeof(Zb), typeof(Xb))  # → false (not defined yet)
```

### Apply Transform Explicitly

```julia
transform(up_x, typeof(Zb))  # → (|↑⟩ + |↓⟩)/√2 as SumKet
```

### Clear All Transforms

```julia
clear_transforms!()  # Remove all registered transforms
```

## Bidirectional Transforms

Transforms are **directional**. If you need to transform in both directions, register both:

```julia
# X → Z (spin-x to spin-z)
define_transform!(Xb, Zb) do idx
    idx == :↑ ? (up_z + down_z)/√2 : (up_z - down_z)/√2
end

# Z → X (spin-z to spin-x)  
define_transform!(Zb, Xb) do idx
    idx == :↑ ? (up_x + down_x)/√2 : (up_x - down_x)/√2
end
```

## Example: Position and Momentum Bases

Here's how you might set up position-momentum transforms (conceptually):

```julia
H = HilbertSpace(:particle)
Hb = Basis(H, :default)

Xb = Basis(H, :x)  # position basis
Pb = Basis(H, :p)  # momentum basis

# In reality, these are continuous bases with Fourier transform relationship
# |p⟩ = ∫ e^{ipx/ℏ}|x⟩ dx / √(2πℏ)
```

!!! note "Discrete vs Continuous"
    QSymbolic.jl currently works with symbolic/discrete states. Continuous bases like position/momentum can be modeled discretely or symbolically, but true continuous integrals are beyond the current scope.
