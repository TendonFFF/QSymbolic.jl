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

For operators with complex or infinite-dimensional action (like ladder operators in Fock space), use `FunctionOperator`:

```julia
F = FockSpace(:mode)
Fb = Basis(F, :n)

# Annihilation operator: â|n⟩ = √n |n-1⟩
â = FunctionOperator(:â, Fb) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end

# Apply to Fock states
n3 = Ket(Fb, 3)
â * n3  # → √3 |2⟩

n0 = Ket(Fb, 0)
â * n0  # → 0 (vacuum annihilated)
```

### When to Use FunctionOperator

Use `FunctionOperator` when:
- The action depends on the state label (like Fock space number)
- The space is infinite-dimensional
- The transformation is more naturally expressed procedurally

Use outer product `Operator` when:
- Building operators from known matrix elements
- Working with finite-dimensional systems
- Combining projectors and ladder operators

## Identity Operator

```julia
# Create identity on a basis
I = IdentityOp(Zb)

I * up    # → |↑⟩
I * down  # → |↓⟩
```

## Basis Dependence

!!! important "Operators Are Basis-Dependent"
    All operators in QSymbolic.jl are associated with a specific basis. This reflects the physical fact that an operator's matrix representation depends on the chosen basis.

```julia
Zb = Basis(H, :z)
Xb = Basis(H, :x)

up_z = Ket(Zb, :↑)
up_x = Ket(Xb, :↑)

# Operator in z-basis
P_z = up_z * up_z'

# Applying to z-basis ket works
P_z * up_z  # → |↑⟩

# Applying to x-basis ket returns symbolic (cross-basis)
P_z * up_x  # → OpKet (symbolic, needs basis transform)
```

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
