# Symbolic Scalars

The symbolic scalar system is powered by [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), providing powerful automatic algebraic simplification and symbolic computation capabilities.

## Overview

| Type | Description |
|:-----|:------------|
| `AbstractSymbolic` | Type alias for `Union{Symbolics.Num, Complex{Symbolics.Num}}` |
| `Sym` | Function to create symbolic variables (returns `Symbolics.Num`) |
| `SymNum` | Wrapper for numeric values (identity function for compatibility) |
| `KroneckerDelta` | Symbolic Kronecker delta δᵢⱼ |

## Types

### Abstract Type

```@docs
AbstractSymbolic
```

### Symbolic Variable

```@docs
Sym
```

### Wrapped Numeric

```@docs
SymNum
```

### Kronecker Delta

```@docs
KroneckerDelta
```

## Evaluation Functions

### Substitution

```@docs
substitute
```

### Numeric Evaluation

```@docs
evaluate
```

### Simplification

```@docs
simplify
```

## Introspection

### Symbol Extraction

```@docs
symbols
```

### Numeric Check

```@docs
is_numeric
```

### Type Assumptions

```@docs
is_real
is_positive
is_negative
is_nonnegative
is_integer
assumptions
```

## Supported Operations

Thanks to Symbolics.jl, the following operations are supported with **automatic algebraic simplification**:

### Arithmetic
- `+`, `-` (binary and unary)
- `*`, `/`, `//`
- `^`

### Functions
- `sqrt`, `√`
- `conj`, `adjoint` (')
- `abs`, `abs2`
- `sin`, `cos`, `tan`
- `exp`, `log`

## Automatic Simplification

One of the key benefits of the Symbolics.jl backend is automatic simplification:

```julia
using QSymbolic

n = Sym(:n, integer=true)

# These simplify automatically:
sqrt(n) * sqrt(n)  # → n
(n - 1) + 1        # → n
n * 1              # → n
n + 0              # → n
```

## Examples

### Basic Usage

```julia
using QSymbolic

# Create symbolic variables
n = Sym(:n, integer=true)
θ = Sym(:θ, real=true)

# Build expressions (auto-simplified by Symbolics.jl)
expr = √n + 1

# Substitute values
result = substitute(expr, :n => 4)
```

### Introspection

```julia
a, b = Sym(:a), Sym(:b)
expr = a^2 + 2*a*b + b^2

symbols(expr)        # → [:a, :b]
is_numeric(expr)     # → false
```

### Using Symbolics.jl Directly

You can also use Symbolics.jl macros directly:

```julia
using QSymbolic

# Using @variables macro (re-exported from Symbolics.jl)
@variables x y::Integer

expr = x^2 + y
```

## See Also

- [Symbolics.jl Documentation](https://docs.sciml.ai/Symbolics/stable/)
- [Guide: Symbolic Scalars](../guide/symbolic.md)
