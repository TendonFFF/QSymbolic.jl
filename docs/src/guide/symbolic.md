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

Symbolic scalars integrate fully with the quantum state system.

### Symbolic Indices

Ket and Bra accept symbolic indices, enabling general representations like |n⟩ for arbitrary n:

```julia
using QSymbolic

F = FockSpace(:fock)
Fb = Basis(F, :n)

# Create symbolic index
n = Sym(:n)

# Ket with symbolic index |n⟩
ket_n = Ket(Fb, n)
ket_n  # → |n⟩

# Bra with symbolic index ⟨n|
bra_n = Bra(Fb, n)
bra_n  # → ⟨n|

# Adjoint preserves symbolic index
ket_n'  # → ⟨n|
```

### Symbolic Weights

Weighted kets and sum states can have symbolic coefficients:

```julia
using QSymbolic

H = HilbertSpace(:H, 2)
Hb = Basis(H, :default)
ψ = Ket(Hb, :ψ)
ϕ = Ket(Hb, :ϕ)

# Symbolic amplitudes
α = Sym(:α)
β = Sym(:β)

# Weighted ket with symbolic coefficient
weighted = α * ψ
weighted  # → (α)|ψ⟩

# Superposition with symbolic amplitudes
state = α * ψ + β * ϕ
state  # → (α)|ψ⟩ + (β)|ϕ⟩

# Extract symbols from weights
symbols(state.weights[1])  # → Set([:α])
```

### Concrete Substitution

You can substitute symbolic values with concrete numbers:

```julia
# Create state with symbolic coefficients
α = Sym(:α)
β = Sym(:β)
state = α * ψ + β * ϕ

# Substitute to get concrete values
concrete_weights = substitute.(state.weights, Ref(Dict(:α => 1/√2, :β => 1/√2)))

# Check all are numeric
all(is_numeric, concrete_weights)  # → true
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
├── Sym         # Symbolic variable (with optional assumptions)
├── SymNum{T}   # Wrapped numeric value
├── SymExpr     # Expression tree (op + args)
└── KroneckerDelta  # δᵢⱼ (1 if i==j, 0 otherwise)
```

All symbolic types are subtypes of `Number`, allowing them to participate in numeric operations.

## Type Assumptions for Sym

Symbolic variables can carry type assumptions that constrain their mathematical properties:

```julia
using QSymbolic

# Create Sym with assumptions
n = Sym(:n, :integer)           # n is an integer
θ = Sym(:θ, :real)              # θ is real
p = Sym(:p, :positive)          # p > 0
E = Sym(:E, :nonnegative)       # E ≥ 0
q = Sym(:q, :negative)          # q < 0
z = Sym(:z, :complex)           # z is complex
```

### Querying Assumptions

```julia
n = Sym(:n, :integer)
p = Sym(:p, :positive)

is_integer(n)       # → true
is_integer(p)       # → false

is_positive(p)      # → true
is_real(p)          # → true (positive implies real)
is_nonnegative(p)   # → true (positive implies nonnegative)

assumptions(n)      # → Set([:integer])
assumptions(p)      # → Set([:positive])
```

### Implicit Assumptions

Some assumptions imply others:
- `:positive` implies `:real` and `:nonnegative`
- `:negative` implies `:real`
- `:integer` is independent (can be positive, negative, or zero)

## KroneckerDelta

The `KroneckerDelta` type represents the Kronecker delta function δᵢⱼ:

```julia
using QSymbolic

n = Sym(:n)
m = Sym(:m)

# Create Kronecker delta
δ = KroneckerDelta(n, m)
δ  # → δ(n,m)

# Evaluates when indices are equal
δ_nn = KroneckerDelta(n, n)
δ_nn  # → 1 (automatically simplifies)

# With concrete values
KroneckerDelta(1, 1)  # → 1
KroneckerDelta(1, 2)  # → 0

# With literal Symbols
KroneckerDelta(:a, :a) |> simplify  # → 1
KroneckerDelta(:a, :b) |> simplify  # → 0

# With symbolic variables (could be equal after substitution)
a = Sym(:a)
b = Sym(:b)
KroneckerDelta(a, b) |> simplify  # → δ(a,b) (not simplified - could be equal)
KroneckerDelta(a, a) |> simplify  # → 1 (same variable)
```

Kronecker deltas arise naturally from inner products of basis states with symbolic indices:

```julia
F = FockSpace(:fock)
Fb = Basis(F, :n)

n = Sym(:n)
m = Sym(:m)

ket_n = Ket(Fb, n)
bra_m = Bra(Fb, m)

# Inner product gives Kronecker delta
bra_m * ket_n  # → δ(m,n)
```

## Custom Contraction Rules

By default, QSymbolic uses orthonormal inner products: `⟨i|j⟩ = δᵢⱼ`. You can define custom contraction rules for specific bases using `define_contraction_rule!`.

### Basic Usage

```julia
using QSymbolic

# Create a custom basis
H = HilbertSpace(:H, 2)
Cb = Basis(H, :coherent)

# Define a custom contraction rule
define_contraction_rule!(typeof(Cb)) do bra_idx, ket_idx
    # Example: coherent states have overlap exp(-|α-β|²/2)
    # For simplicity, returning symbolic expression
    if bra_idx isa AbstractSymbolic || ket_idx isa AbstractSymbolic
        # Keep symbolic
        return KroneckerDelta(bra_idx, ket_idx)
    else
        # Compute numerically
        return bra_idx == ket_idx ? 1 : 0
    end
end
```

### Working with Symbolic Indices

!!! warning "Important: Use `isequal()` instead of `==`"
    When writing contraction rules that might receive symbolic values (`Symbolics.Num`),
    **never use `==` for comparisons**. The `==` operator on symbolic values returns
    a symbolic expression, not a boolean, which will cause errors in `if` statements.

    Instead, use:
    - `isequal(a, b)` - Returns `true` if `a` and `b` are identical objects
    - `a isa AbstractSymbolic` - Check if a value is symbolic
    - `a isa Symbol` - Check if a value is a literal Julia Symbol

### Example: Jaynes-Cummings Model

Here's a complete example for a composite basis in a Jaynes-Cummings system:

```julia
using QSymbolic

# Create spaces
S_cavity, B_cavity = FockSpace(:cavity)
S_dot = HilbertSpace(:dot, 3)

# Create composite basis
S_system = S_cavity ⊗ S_dot
B_JC = Basis(S_system, :JC)

# Define contraction rule for JC basis
define_contraction_rule!(typeof(B_JC)) do bra_ind, ket_ind
    m, ξ = bra_ind  # (photon number, dot state)
    n, σ = ket_ind

    # Handle special ground state symbol
    # Use isequal() NOT == for symbol comparison!
    if isequal(m, :g)
        if isequal(n, :g)
            return 1
        elseif n isa Symbol
            return 0  # Different literal symbols
        end
        return KroneckerDelta(m, n) |> simplify
    end

    # Check if values are symbolic
    m_symbolic = m isa AbstractSymbolic
    n_symbolic = n isa AbstractSymbolic
    ξ_symbolic = ξ isa AbstractSymbolic
    σ_symbolic = σ isa AbstractSymbolic

    # Photon number overlap
    if !m_symbolic && !n_symbolic
        num = isequal(m, n) ? 1 : 0
    else
        num = KroneckerDelta(m, n) |> simplify
    end

    # Dot state overlap
    if !ξ_symbolic && !σ_symbolic
        par = isequal(ξ, σ) ? 1 : 0
    else
        par = KroneckerDelta(ξ, σ) |> simplify
    end

    return num * par
end

# Now inner products work correctly
g0 = Ket(B_JC, (Sym(:n), Sym(:σ)))
g0' * g0  # → 1 (same symbolic indices)

k1 = Ket(B_JC, (1, 2))
k2 = Ket(B_JC, (1, 3))
k1' * k2  # → 0 (different second index)
```

### Contraction Rule API

| Function | Description |
|:---------|:------------|
| `define_contraction_rule!(typeof(basis), rule)` | Define custom rule for a basis type |
| `has_contraction_rule(typeof(basis))` | Check if custom rule exists |
| `apply_contraction_rule(typeof(basis), i, j)` | Manually apply rule |
| `clear_contraction_rules!()` | Remove all custom rules |

## API Summary

| Function | Description |
|:---------|:------------|
| `Sym(:name)` | Create symbolic variable |
| `Sym(:name, :assumption)` | Create symbolic variable with assumption |
| `SymNum(value)` | Wrap numeric value |
| `KroneckerDelta(i, j)` | Create Kronecker delta δᵢⱼ |
| `substitute(expr, pairs...)` | Replace symbols with values |
| `evaluate(expr)` | Compute numeric result |
| `simplify(expr)` | Apply algebraic simplifications |
| `symbols(expr)` | Get set of symbol names |
| `is_numeric(expr)` | Check if fully numeric |
| `is_real(sym)` | Check if Sym has real assumption |
| `is_positive(sym)` | Check if Sym has positive assumption |
| `is_integer(sym)` | Check if Sym has integer assumption |
| `assumptions(sym)` | Get set of assumptions |
