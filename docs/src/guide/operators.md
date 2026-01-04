# Operators

QSymbolic.jl provides a comprehensive system for quantum operators, including outer product operators, operator algebra, and function-defined operators.

## Outer Product Operators

The most natural way to build quantum operators is via outer products `|ψ⟩⟨ϕ|`. When applied to a state:

```math
(|ψ⟩⟨ϕ|)|χ⟩ = ⟨ϕ|χ⟩ |ψ⟩
```

```julia
using QSymbolic

H = HilbertSpace(:spin, 2)
Zb = Basis(H, :z)

up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Create operator via outer product
P_up = up * up'      # |↑⟩⟨↑|

# Display shows Dirac notation
P_up  # → |↑⟩⟨↑|

# Apply to states
P_up * up    # → |↑⟩ (eigenstate)
P_up * down  # → 0   (orthogonal)
```

## Projectors

A **projector** onto state `|ψ⟩` is `P = |ψ⟩⟨ψ|`. It satisfies:
- `P|ψ⟩ = |ψ⟩` (leaves the state unchanged)
- `P|ϕ⟩ = 0` for orthogonal states
- `P² = P` (idempotent)

```julia
# Projector onto |↑⟩
P_up = up * up'

P_up * up                # → |↑⟩
P_up * down              # → 0
(P_up * P_up) * up       # → |↑⟩ (P² = P)
```

## Ladder (Raising/Lowering) Operators

Operators like `|↑⟩⟨↓|` raise or lower states:

```julia
# Raising operator σ₊ = |↑⟩⟨↓|
σ_plus = up * down'

σ_plus * down  # → |↑⟩ (raises |↓⟩ to |↑⟩)
σ_plus * up    # → 0   (can't raise |↑⟩ further)

# Lowering operator σ₋ = |↓⟩⟨↑| = (σ₊)†
σ_minus = σ_plus'

σ_minus * up    # → |↓⟩
σ_minus * down  # → 0
```

## Operator Algebra

Operators can be combined using standard arithmetic:

### Addition and Subtraction

```julia
P_up = up * up'
P_down = down * down'

# Identity operator: 𝕀 = |↑⟩⟨↑| + |↓⟩⟨↓|
I = P_up + P_down

I * up    # → |↑⟩
I * down  # → |↓⟩

# Pauli Z: σz = |↑⟩⟨↑| - |↓⟩⟨↓|
σz = P_up - P_down

σz * up    # → |↑⟩  (eigenvalue +1)
σz * down  # → -|↓⟩ (eigenvalue -1)
```

### Scalar Multiplication

```julia
# Scale an operator
half_σz = (1/2) * σz

half_σz * up  # → 0.5|↑⟩
```

### Operator Products

```julia
# Product of operators: (ÂB̂)|ψ⟩ = Â(B̂|ψ⟩)
σ_plus = up * down'
σ_minus = down * up'

# σ₊σ₋ = |↑⟩⟨↓|↓⟩⟨↑| = |↑⟩⟨↑|
product = σ_plus * σ_minus

product * up    # → |↑⟩
product * down  # → 0
```

## Adjoint (Hermitian Conjugate)

The adjoint of `|ψ⟩⟨ϕ|` is `|ϕ⟩⟨ψ|`:

```julia
A = up * down'   # |↑⟩⟨↓|
A'               # → |↓⟩⟨↑|
```

For sum operators, the adjoint distributes:
```julia
σz = P_up - P_down
σz'  # → |↑⟩⟨↑| - |↓⟩⟨↓| (σz is Hermitian)
```

## Applying Operators to Superpositions

Operators act linearly on superpositions:

```julia
# Superposition |ψ⟩ = (|↑⟩ + |↓⟩)/√2
ψ = (up + down) / √2

# σz|ψ⟩ = (|↑⟩ - |↓⟩)/√2
σz * ψ  # → 0.707|↑⟩ - 0.707|↓⟩
```

## Function-Based Operators

For operators with complex or infinite-dimensional action (like ladder operators in Fock space), use `FunctionOperator`. This is particularly useful when:

- The action depends on the state label (like Fock space number)
- The space is infinite-dimensional
- The transformation is more naturally expressed procedurally

### Basic Syntax

```julia
F, Fb = FockSpace(:mode)

# Define action function
annihilate(ket::Ket{B}) where B = √(ket.index) * Ket{B}(ket.index - 1)

# Create function operator
â = FunctionOperator(annihilate, Fb, name=:â)

# Apply to Fock states
n3 = Ket(Fb, 3)
â * n3  # → √3 |2⟩

n0 = Ket(Fb, 0)
â * n0  # → 0 (vacuum annihilated)
```

### With Do-Block Syntax

```julia
F, Fb = FockSpace(:mode)

# Using do-block
â = FunctionOperator(Fb, name=:â) do ket
    n = ket.index
    n == 0 ? 0 : √n * Ket{typeof(Fb)}(n - 1)
end
```

### With Adjoint Action

For operators where you also need the adjoint (like creation operator for annihilation):

```julia
F, Fb = FockSpace(:mode)

# Define both actions
annihilate(ket::Ket{B}) where B = √(ket.index) * Ket{B}(ket.index - 1)
create(ket::Ket{B}) where B = √(ket.index + 1) * Ket{B}(ket.index + 1)

# Create with adjoint
â = FunctionOperator(annihilate, Fb, adjoint_action=create, name=:â)

# Now adjoint works
â'  # Creation operator â†

# Verify â†â|n⟩ = n|n⟩ (number operator)
n = Sym(:n, :nonnegative, :integer)
ket_n = Ket(Fb, n)
â' * (â * ket_n)  # → n|n⟩
```

### Working with Symbolic Indices

Function operators work seamlessly with symbolic indices:

```julia
F, Fb = FockSpace(:mode)

# Action handles symbolic indices
annihilate(ket::Ket{B}) where B = √(ket.index) * Ket{B}(ket.index - 1)
â = FunctionOperator(annihilate, Fb, name=:â)

# Symbolic index
d = Sym(:d, :nonnegative, :integer)
ket_d = Ket(Fb, d)

â * ket_d  # → √d |d-1⟩
```

### Operators on Composite States

Function operators automatically handle tensor products:

```julia
S_cavity, B_cavity = FockSpace(:cavity)
S_dot, B_dot = HilbertSpace(:dot, 3)

# Annihilation operator on cavity
annihilate(ket::Ket{B}) where B = √(ket.index) * Ket{B}(ket.index - 1)
a = FunctionOperator(annihilate, B_cavity, name=:a)

# Apply to product state |d⟩ ⊗ |g⟩
d = Sym(:d)
product_state = Ket(B_cavity, d) ⊗ Ket(B_dot, :g)

a * product_state  # → √d |d-1⟩ ⊗ |g⟩
```

## Identity Operator

```julia
H = HilbertSpace(:H, 2)

# Create identity on a space
I = Identity(H)

H, Hb = HilbertSpace(:H, 2)
ψ = Ket(Hb, :ψ)
I * ψ  # → |ψ⟩
```

## Operator Types Summary

| Type | Description | Example |
|:-----|:------------|:--------|
| `Outer` | Single outer product `\|ψ⟩⟨ϕ\|` | `up * down'` |
| `Operator` | Sum of weighted outers | `σz = P_up - P_down` |
| `Identity` | Identity on space | `Identity(H)` |
| `FunctionOperator` | User-defined action | `FunctionOperator(action, basis)` |

## Complete Example: Pauli Matrices

```julia
using QSymbolic

H = HilbertSpace(:spin, 2)
Zb = Basis(H, :z)

up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Build Pauli matrices from outer products
σx = up * down' + down * up'         # |↑⟩⟨↓| + |↓⟩⟨↑|
σy = -im * up * down' + im * down * up'  # -i|↑⟩⟨↓| + i|↓⟩⟨↑|
σz = up * up' - down * down'         # |↑⟩⟨↑| - |↓⟩⟨↓|

# Verify eigenvalues
σz * up    # → |↑⟩   (eigenvalue +1)
σz * down  # → -|↓⟩  (eigenvalue -1)

σx * up    # → |↓⟩
σx * down  # → |↑⟩

# Apply to superposition
ψ = (up + down) / √2
σz * ψ     # → (|↑⟩ - |↓⟩)/√2
```

## Complete Example: Fock Space Ladder Operators

```julia
using QSymbolic

F, Fb = FockSpace(:mode)

# Define annihilation and creation
annihilate(ket::Ket{B}) where B = √(ket.index) * Ket{B}(ket.index - 1)
create(ket::Ket{B}) where B = √(ket.index + 1) * Ket{B}(ket.index + 1)

â = FunctionOperator(annihilate, Fb, adjoint_action=create, name=:â)
â† = â'

# Number operator N̂ = â†â
# Apply to |n⟩
n = Sym(:n, :nonnegative, :integer)
ket_n = Ket(Fb, n)

â * ket_n       # → √n |n-1⟩
â† * ket_n      # → √(n+1) |n+1⟩
â† * (â * ket_n)  # → n|n⟩ (eigenvalue equation)
```

## See Also

- [Getting Started](@ref) - Basic quantum state operations
- [Composite Systems](@ref) - Operators on tensor products
- [Basis Transforms](@ref) - Cross-basis operator application
