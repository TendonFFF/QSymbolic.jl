# Transforms

Basis transformations allow you to express states in different representations. QSymbolic.jl maintains a registry of user-defined transforms and applies them automatically when computing cross-basis inner products.

## Overview

The transform system works as follows:
1. **Register** transforms between bases using `define_transform!`
2. **Query** whether transforms exist with `has_transform`
3. **Apply** transforms explicitly with `transform` or let inner products use them automatically

For composite systems, transforms **factorize automatically** — if you define transforms for each subsystem, the composite transform is derived without additional registration.

## Registering Transforms

```@docs
define_transform!
```

## Querying Transforms

```@docs
has_transform
```

## Applying Transforms

```@docs
transform
```

## Utilities

```@docs
clear_transforms!
```

## Examples

### Basic Transform

```julia
using QSymbolic

H = HilbertSpace(:spin, 2)
Hb = Basis(H, :default)
Zb = Basis(H, :z)
Xb = Basis(H, :x)

up_z = Ket(Zb, :↑)
down_z = Ket(Zb, :↓)

# Register X → Z transform
define_transform!(Xb, Zb) do idx
    idx == :↑ ? (up_z + down_z)/√2 : (up_z - down_z)/√2
end

# Check availability
has_transform(typeof(Xb), typeof(Zb))  # true

# Apply transform
up_x = Ket(Xb, :↑)
transform(up_x, typeof(Zb))  # (|↑⟩ + |↓⟩)/√2

# Inner products now work automatically
up_z' * up_x  # 1/√2
```

### Composite Transform (Factorized)

```julia
H_A, H_B = HilbertSpace(:A, 2), HilbertSpace(:B, 2)
Za, Xa = Basis(H_A, :z), Basis(H_A, :x)
Zb, Xb = Basis(H_B, :z), Basis(H_B, :x)

# Define individual transforms (same pattern as above)
define_transform!(Xa, Za) do idx
    # ... x-to-z for subsystem A
end
define_transform!(Xb, Zb) do idx
    # ... x-to-z for subsystem B  
end

# Composite transform is automatic!
has_transform(typeof(Xa ⊗ Xb), typeof(Za ⊗ Zb))  # true
```
