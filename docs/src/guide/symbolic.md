# Symbolic Scalars

QSymbolic.jl includes a symbolic scalar system for lazy arithmetic evaluation. This enables truly symbolic quantum computations where numerical values are substituted only when needed.

## Why Symbolic Scalars?

In many quantum mechanical calculations, you want to keep expressions symbolic:
- Derive general formulas with parameters like `n`, `θ`, `ω`
- Substitute specific values later
- Avoid floating-point errors in intermediate steps
- Generate readable analytical expressions

## Creating Symbolic Variables

Use `Sym(:name)` to create a symbolic variable:

```julia
using QSymbolic

n = Sym(:n)
θ = Sym(:θ)
ω = Sym(:ω)

# Display shows the symbol name
n  # → n
```

## Building Expressions

Arithmetic operations on symbolic variables build expression trees (no evaluation happens):

```julia
n = Sym(:n)

# Basic arithmetic
n + 1        # → (n + 1)
n - 2        # → (n - 2)
n * 3        # → n·3
n / 4        # → (n/4)
n^2          # → n^2

# Square root
√n           # → √(n)
sqrt(n)      # → √(n)

# Combinations
√n + 1       # → (√(n) + 1)
n^2 + 2*n    # → ((n^2) + (2·n))
```

## Supported Operations

| Operation | Syntax | Example |
|:----------|:-------|:--------|
| Addition | `+` | `n + 1` |
| Subtraction | `-` | `n - 1` |
| Negation | `-` | `-n` |
| Multiplication | `*` | `n * 2` |
| Division | `/`, `//` | `n / 2` |
| Power | `^` | `n^2` |
| Square root | `√`, `sqrt` | `√n` |
| Conjugate | `conj`, `'` | `conj(z)` |
| Absolute value | `abs`, `abs2` | `abs(z)` |
| Trigonometric | `sin`, `cos`, `tan` | `sin(θ)` |
| Exponential | `exp`, `log` | `exp(x)` |

## Substitution

Replace symbolic variables with values using `substitute`:

```julia
n = Sym(:n)
expr = √n + 1

# Substitute n → 4
result = substitute(expr, :n => 4)
result  # → (√(4) + 1)  (still symbolic, but with numeric value)
```

Multiple substitutions:
```julia
a = Sym(:a)
b = Sym(:b)
expr = a^2 + b^2

substitute(expr, :a => 3, :b => 4)  # → ((3^2) + (4^2))
```

Partial substitution:
```julia
expr = a * b + a
substitute(expr, :a => 2)  # → ((2·b) + 2)  (b remains symbolic)
```

## Evaluation

Convert a fully-substituted expression to a numeric value with `evaluate`:

```julia
n = Sym(:n)
expr = √n + 1

# Substitute then evaluate
result = substitute(expr, :n => 4)
evaluate(result)  # → 3.0

# Chained
substitute(expr, :n => 9) |> evaluate  # → 4.0
```

!!! warning "All Symbols Must Be Substituted"
    `evaluate` throws an error if any symbolic variables remain:
    ```julia
    evaluate(Sym(:n))  # Error: Cannot evaluate: unsubstituted symbol n
    ```

## Introspection

### List Symbols in Expression

```julia
n = Sym(:n)
m = Sym(:m)
expr = n^2 + 2*n*m + m^2

symbols(expr)  # → Set([:n, :m])
```

### Check if Numeric

```julia
n = Sym(:n)
expr = √n + 1

is_numeric(n)                    # → false
is_numeric(SymNum(5))            # → true
is_numeric(expr)                 # → false (contains n)
is_numeric(substitute(expr, :n => 4))  # → true (all substituted)
```

## Simplification

Apply basic algebraic simplifications with `simplify`:

```julia
n = Sym(:n)

# Identity simplifications
simplify(n * 1)   # → n
simplify(n * 0)   # → 0
simplify(n + 0)   # → n
simplify(n^1)     # → n
simplify(n^0)     # → 1

# Numeric folding
simplify(SymNum(2) + SymNum(3))  # → 5
simplify(SymNum(2) * SymNum(3))  # → 6
```

## Use with Quantum States

Symbolic scalars integrate with the quantum state system:

```julia
using QSymbolic

H = HilbertSpace(:H, 2)
Zb = Basis(H, :z)
up = BasisKet(Zb, :↑)
down = BasisKet(Zb, :↓)

# Symbolic coefficient
θ = Sym(:θ)

# State with symbolic amplitude
# (In general form - for actual use you'd define the arithmetic)
coeff = cos(θ)

# Build superposition with symbolic weights
# Note: Full integration depends on state arithmetic supporting AbstractSymbolic
```

## Example: Parametric Expressions

```julia
using QSymbolic

# Define parameters
n = Sym(:n)
ω = Sym(:ω)
t = Sym(:t)

# Build expression for harmonic oscillator energy
E_n = ω * (n + 1//2)
E_n  # → ω·(n + 1//2)

# Substitute specific values
E_0 = substitute(E_n, :n => 0)
E_0  # → ω·(0 + 1//2)

E_0_numeric = substitute(E_0, :ω => 2π)
evaluate(E_0_numeric)  # → π ≈ 3.14159...
```

## Type Hierarchy

```
AbstractSymbolic <: Number
├── Sym         # Symbolic variable
├── SymNum{T}   # Wrapped numeric value
└── SymExpr     # Expression tree (op + args)
```

All symbolic types are subtypes of `Number`, allowing them to participate in numeric operations.

## API Summary

| Function | Description |
|:---------|:------------|
| `Sym(:name)` | Create symbolic variable |
| `SymNum(value)` | Wrap numeric value |
| `substitute(expr, pairs...)` | Replace symbols with values |
| `evaluate(expr)` | Compute numeric result |
| `simplify(expr)` | Apply algebraic simplifications |
| `symbols(expr)` | Get set of symbol names |
| `is_numeric(expr)` | Check if fully numeric |
