# Symbolic Scalars

The symbolic scalar system provides lazy arithmetic evaluation, enabling truly symbolic quantum computations.

## Overview

| Type | Description |
|:-----|:------------|
| `AbstractSymbolic` | Abstract supertype for all symbolic scalars |
| `Sym` | Symbolic variable |
| `SymNum` | Wrapped numeric value |
| `SymExpr` | Expression tree |

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

### Expression Tree

```@docs
SymExpr
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

## Supported Operations

The following operations are defined for `AbstractSymbolic` types and build expression trees:

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

## Examples

### Basic Usage

```julia
using QSymbolic

# Create symbolic variables
n = Sym(:n)
θ = Sym(:θ)

# Build expressions
expr = √n + 1
trig = sin(θ)^2 + cos(θ)^2

# Substitute and evaluate
result = substitute(expr, :n => 4)
evaluate(result)  # → 3.0
```

### Introspection

```julia
a, b = Sym(:a), Sym(:b)
expr = a^2 + 2*a*b + b^2

symbols(expr)        # → Set([:a, :b])
is_numeric(expr)     # → false
is_numeric(substitute(expr, :a => 1, :b => 2))  # → true
```

### Simplification

```julia
n = Sym(:n)

simplify(n * 1)  # → n
simplify(n * 0)  # → 0
simplify(n + 0)  # → n
simplify(n^0)    # → 1
```
