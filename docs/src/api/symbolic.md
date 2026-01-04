# Symbolic Scalars

QSymbolic.jl uses [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) as its backend for symbolic computation.

## Overview

| Type/Function | Description |
|:--------------|:------------|
| `AbstractSymbolic` | Type alias for Symbolics.jl types |
| `Sym` | Create symbolic variable with assumptions |
| `@variables`, `@syms` | Re-exported from Symbolics.jl |
| `KroneckerDelta` | Kronecker delta δᵢⱼ |
| `substitute` | Replace symbols with values |
| `evaluate` | Compute numeric result |
| `simplify` | Algebraic simplification |

## Types

### AbstractSymbolic

Type alias for symbolic types:

```julia
const AbstractSymbolic = Union{Symbolics.Num, Complex{Symbolics.Num}}
```

### Symbolic Variable

```@docs
Sym
```

### Kronecker Delta

```@docs
KroneckerDelta
```

## Creating Variables

```julia
# QSymbolic API
n = Sym(:n)                        # Default (real)
m = Sym(:m, :integer)              # With assumption
θ = Sym(:θ, :real, :positive)      # Multiple assumptions

# Symbolics.jl macros (re-exported)
@variables x y z
@syms a b c
@variables n::Integer
```

## Assumptions

Available assumptions for `Sym`:
- `:real` - real number
- `:positive` - strictly positive (implies real)
- `:negative` - strictly negative (implies real)
- `:nonnegative` - non-negative (implies real)
- `:integer` - integer value
- `:complex` - complex number

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

```@docs
symbols
is_numeric
is_real
is_positive
is_integer
assumptions
```

## KroneckerDelta Behavior

```julia
# Concrete values
KroneckerDelta(1, 1)  # → 1
KroneckerDelta(1, 2)  # → 0

# Literal Symbols
KroneckerDelta(:a, :a) |> simplify  # → 1
KroneckerDelta(:a, :b) |> simplify  # → 0

# Symbolic variables (may be equal after substitution)
a, b = Sym(:a), Sym(:b)
KroneckerDelta(a, b) |> simplify    # → δ(a,b) (stays symbolic)
KroneckerDelta(a, a) |> simplify    # → 1 (same variable)
```

## Examples

```julia
using QSymbolic

# Create symbolic variables
n = Sym(:n, :nonnegative, :integer)
ω = Sym(:ω, :positive)

# Build expressions
E_n = ω * (n + 1//2)  # Harmonic oscillator energy

# Substitute values
E_0 = substitute(E_n, :n => 0)  # → ω/2
E_numeric = substitute(E_0, :ω => 2π) |> evaluate  # → π

# Use in quantum states
F, Fb = FockSpace(:mode)
ket_n = Ket(Fb, n)
ket_m = Ket(Fb, Sym(:m))
ket_m' * ket_n  # → δ(m,n)
```
