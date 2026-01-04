# Symbolic Scalars

QSymbolic.jl uses [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) as its backend for symbolic computation. A thin compatibility layer provides a familiar API while leveraging Symbolics.jl's powerful computer algebra system.

## Why Symbolic Scalars?

In quantum mechanical calculations, you often want to keep expressions symbolic:
- Derive general formulas with parameters like `n`, `θ`, `ω`
- Substitute specific values later
- Avoid floating-point errors in intermediate steps
- Generate readable analytical expressions

## Creating Symbolic Variables

### Using `Sym` (QSymbolic API)

The `Sym` constructor creates symbolic variables with optional type assumptions:

```julia
using QSymbolic

n = Sym(:n)           # Generic symbolic variable
θ = Sym(:θ)           # Greek letters work too
ω = Sym(:ω, :real)    # Real-valued

# Display shows the symbol name
n  # → n
```

### Using Symbolics.jl Macros (Re-exported)

QSymbolic re-exports `@variables` and `@syms` from Symbolics.jl:

```julia
using QSymbolic

# Create multiple variables at once
@variables x y z
@syms a b c

# With type annotations
@variables t::Real
@variables n::Integer
```

## Type Assumptions

Symbolic variables can carry assumptions that affect simplification and conjugation:

```julia
# Positional syntax
n = Sym(:n, :integer)           # n is an integer
θ = Sym(:θ, :real)              # θ is real
p = Sym(:p, :positive)          # p > 0 (implies real)
E = Sym(:E, :nonnegative)       # E ≥ 0

# Keyword syntax
k = Sym(:k, integer=true)
r = Sym(:r, real=true, positive=true)

# Multiple assumptions
m = Sym(:m, :real, :positive, :integer)
```

### Querying Assumptions

```julia
n = Sym(:n, :integer)
p = Sym(:p, :positive)

is_integer(n)       # → true
is_positive(p)      # → true
is_real(p)          # → true (positive implies real)
is_nonnegative(p)   # → true (positive implies nonnegative)

assumptions(n)      # → Set([:integer])
```

### Effect on Conjugation

Assumptions affect how the adjoint (complex conjugate) behaves:

```julia
z = Sym(:z)                    # Complex by default
z'                             # → conj(z)

r = Sym(:r, :real)
r'                             # → r (real variables self-conjugate)
```

## Building Expressions

Arithmetic operations build expression trees automatically:

```julia
n = Sym(:n)

# Basic arithmetic
n + 1        # → n + 1
n - 2        # → n - 2
n * 3        # → 3n
n / 4        # → n/4
n^2          # → n²

# Square root
√n           # → √n
sqrt(n)      # → √n

# Combinations
√n + 1       # → 1 + √n
n^2 + 2n     # → n² + 2n
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
| Conjugate | `conj` | `conj(z)` |
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
result  # → 1 + √4 = 3

# Multiple substitutions
@syms a b
expr = a^2 + b^2
substitute(expr, :a => 3, :b => 4)  # → 25

# Partial substitution
expr = a * b + a
substitute(expr, :a => 2)  # → 2b + 2 (b remains symbolic)
```

## Evaluation

Convert expressions to numeric values with `evaluate`:

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
    `evaluate` throws an error if symbolic variables remain:
    ```julia
    evaluate(Sym(:n))  # Error: Cannot evaluate expression with free symbols
    ```

## Simplification

Apply algebraic simplifications with `simplify`:

```julia
n = Sym(:n)

# Identity simplifications
simplify(n * 1)   # → n
simplify(n * 0)   # → 0
simplify(n + 0)   # → n
simplify(n^1)     # → n
simplify(n^0)     # → 1

# Numeric folding
simplify(2 + 3)   # → 5

# Expression simplification
simplify(n + n)   # → 2n
```

## Introspection

### List Symbols in Expression

```julia
@syms n m
expr = n^2 + 2*n*m + m^2

symbols(expr)  # → Set([:n, :m])
```

### Check if Numeric

```julia
n = Sym(:n)
expr = √n + 1

is_numeric(n)                        # → false
is_numeric(5)                        # → true
is_numeric(expr)                     # → false (contains n)
is_numeric(substitute(expr, :n => 4))  # → true
```

## KroneckerDelta

The `KroneckerDelta` type represents the Kronecker delta function δᵢⱼ:

```math
\delta_{ij} = \begin{cases} 1 & \text{if } i = j \\ 0 & \text{if } i \neq j \end{cases}
```

```julia
using QSymbolic

n = Sym(:n)
m = Sym(:m)

# Create Kronecker delta
δ = KroneckerDelta(n, m)
δ  # → δ(n,m)

# Same variable simplifies to 1
KroneckerDelta(n, n) |> simplify  # → 1

# Concrete values evaluate immediately
KroneckerDelta(1, 1)  # → 1
KroneckerDelta(1, 2)  # → 0

# Literal Symbols (not Sym) compare directly
KroneckerDelta(:a, :a) |> simplify  # → 1
KroneckerDelta(:a, :b) |> simplify  # → 0

# Different symbolic variables remain symbolic
a = Sym(:a)
b = Sym(:b)
KroneckerDelta(a, b) |> simplify  # → δ(a,b) (could be equal after substitution)
```

### KroneckerDelta in Inner Products

Kronecker deltas arise naturally from inner products with symbolic indices:

```julia
H, Hb = HilbertSpace(:H, 2)

n = Sym(:n)
m = Sym(:m)

ket_n = Ket(Hb, n)
ket_m = Ket(Hb, m)

# Inner product gives Kronecker delta
ket_m' * ket_n  # → δ(m,n)
ket_n' * ket_n  # → 1 (same index)
```

### Multi-Index Kronecker Deltas

For composite indices, KroneckerDelta products merge:

```julia
n, m = Sym(:n), Sym(:m)
σ, ξ = Sym(:σ), Sym(:ξ)

# Product of deltas
δ1 = KroneckerDelta(n, m)
δ2 = KroneckerDelta(σ, ξ)

δ1 * δ2  # → δ(n,m)·δ(σ,ξ)
```

## Use with Quantum States

### Symbolic Indices

Kets accept symbolic indices for general representations:

```julia
using QSymbolic

S, B = FockSpace(:fock)

# Create symbolic index
n = Sym(:n)

# Ket with symbolic index |n⟩
ket_n = Ket(B, n)
ket_n  # → |n⟩

# Adjoint preserves symbolic index
ket_n'  # → ⟨n|

# Inner products
ket_n' * ket_n  # → 1
m = Sym(:m)
Ket(B, m)' * ket_n  # → δ(m,n)
```

### Symbolic Weights

Weighted kets and superpositions can have symbolic coefficients:

```julia
H, Hb = HilbertSpace(:H, 2)
ψ = Ket(Hb, :ψ)
ϕ = Ket(Hb, :ϕ)

# Symbolic amplitudes
α = Sym(:α)
β = Sym(:β)

# Weighted ket
α * ψ  # → α·|ψ⟩

# Superposition with symbolic amplitudes
state = α * ψ + β * ϕ  # → α·|ψ⟩ + β·|ϕ⟩
```

### Example: Parametric Expressions

```julia
using QSymbolic

# Define parameters
n = Sym(:n, :nonnegative, :integer)
ω = Sym(:ω, :positive)

# Harmonic oscillator energy
E_n = ω * (n + 1//2)
E_n  # → ω*(n + 1/2)

# Substitute specific values
E_0 = substitute(E_n, :n => 0)
E_0  # → ω/2

E_0_numeric = substitute(E_0, :ω => 2π)
evaluate(E_0_numeric)  # → π ≈ 3.14159...
```

## Type System

QSymbolic uses Symbolics.jl types internally:

```julia
# AbstractSymbolic is an alias for Symbolics types
const AbstractSymbolic = Union{Symbolics.Num, Complex{Symbolics.Num}}

# Check if a value is symbolic
x isa AbstractSymbolic

# Sym creates Symbolics.Num variables
n = Sym(:n)
typeof(n)  # → Num
```

### Compatibility Types

For backward compatibility, these constructors are provided:

| Constructor | Behavior |
|:------------|:---------|
| `Sym(:name)` | Creates `Symbolics.Num` variable |
| `SymNum(x)` | Returns `x` unchanged (passthrough) |
| `SymExpr(op, args...)` | Builds expression via arithmetic |

## API Summary

| Function | Description |
|:---------|:------------|
| `Sym(:name)` | Create symbolic variable |
| `Sym(:name, :assumption...)` | Create with assumptions |
| `@variables x y z` | Create multiple variables (Symbolics.jl) |
| `@syms a b c` | Create multiple variables (Symbolics.jl) |
| `KroneckerDelta(i, j)` | Kronecker delta δᵢⱼ |
| `substitute(expr, pairs...)` | Replace symbols with values |
| `evaluate(expr)` | Compute numeric result |
| `simplify(expr)` | Apply algebraic simplifications |
| `symbols(expr)` | Get set of symbol names |
| `is_numeric(expr)` | Check if fully numeric |
| `is_real(sym)` | Check real assumption |
| `is_positive(sym)` | Check positive assumption |
| `is_integer(sym)` | Check integer assumption |
| `assumptions(sym)` | Get set of assumptions |

## See Also

- [Custom Contraction Rules](@ref) - Define non-orthonormal inner products
- [Getting Started](@ref) - Basic quantum state operations
