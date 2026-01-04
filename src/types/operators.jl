# Operator types: ALL operator struct definitions
# Outer, Operator, Identity, FunctionOperator, OperatorSum

# Exports
export Outer, Operator, Identity, FunctionOperator
export OperatorSum
export tr

# space() - get the space an operator acts on
space(::AbstractOperator{S}) where S = S

# ============== Outer Product |ψ⟩⟨ϕ| ==============

@doc """
    Outer{S<:AbstractSpace}(ket, bra)
    ket * bra'  (automatically creates Outer)

Basic outer product operator |ψ⟩⟨ϕ|. Allows cross-basis operators
where the ket and bra can be in different bases, as long as they're
in the same space.

Application uses bra-ket contraction:
    (|ψ⟩⟨ϕ|)|χ⟩ = ⟨ϕ|χ⟩ |ψ⟩

# Examples
```julia
H, Hb = HilbertSpace(:H, 2)
up = Ket(Hb, :↑)
down = Ket(Hb, :↓)

# Projector onto |↑⟩
P_up = up * up'  # |↑⟩⟨↑|

# Apply: P_up|↑⟩ = |↑⟩, P_up|↓⟩ = 0
P_up * up   # → |↑⟩
P_up * down # → 0

# Ladder operator |↑⟩⟨↓|
σ_plus = up * down'
σ_plus * down  # → |↑⟩
```

See also: [`Operator`](@ref), [`Identity`](@ref), [`FunctionOperator`](@ref)
"""
struct Outer{S<:AbstractSpace} <: AbstractOperator{S}
    ket::AbstractKet
    bra::AbstractBra
    
    function Outer(ket::AbstractKet, bra::AbstractBra)
        # Verify same space
        ket_space = space(basis(ket))
        bra_space = space(basis(bra))
        ket_space == bra_space || throw(ArgumentError("Ket and bra must be in the same space"))
        new{ket_space}(ket, bra)
    end
end

# Get space of an Outer
space(::Outer{S}) where S = S

# Create outer product from ket * bra'
Base.:*(ket::AbstractKet, bra::AbstractBra) = Outer(ket, bra)

# Display
function Base.show(io::IO, op::Outer)
    print(io, _ket_str(op.ket), _bra_str(op.bra))
end

# Helper functions for display
_ket_str(k::Ket) = "|$(k.index)⟩"
_ket_str(k::WeightedKet) = "$(k.weight)·|$(k.ket.index)⟩"
_ket_str(k::ProductKet) = "|" * join([string(ki.index) for ki in k.kets], "⊗") * "⟩"
_ket_str(k::SumKet) = "(" * join(["$(w)·|$(ki.index)⟩" for (ki, w) in zip(k.kets, k.weights)], " + ") * ")"
_ket_str(k::AbstractKet) = "|ket⟩"

_bra_str(b::Bra) = "⟨$(b.index)|"
_bra_str(b::WeightedBra) = "$(b.weight)·⟨$(b.bra.index)|"
_bra_str(b::ProductBra) = "⟨" * join([string(bi.index) for bi in b.bras], "⊗") * "|"
_bra_str(b::SumBra) = "(" * join(["$(w)·⟨$(bi.index)|" for (bi, w) in zip(b.bras, b.weights)], " + ") * ")"
_bra_str(b::AbstractBra) = "⟨bra|"

# ============== Operator Container ==============

@doc """
    Operator{S<:AbstractSpace}(outers, weights)

Container type representing a sum of weighted outer products:
    Σᵢ wᵢ |ψᵢ⟩⟨ϕᵢ|

Application decomposes to bra-ket arithmetic:
    (Σᵢ wᵢ |ψᵢ⟩⟨ϕᵢ|)|χ⟩ = Σᵢ wᵢ ⟨ϕᵢ|χ⟩ |ψᵢ⟩

# Examples
```julia
H, Hb = HilbertSpace(:H, 2)
up = Ket(Hb, :↑)
down = Ket(Hb, :↓)

# σ_x = |↑⟩⟨↓| + |↓⟩⟨↑|
σ_x = Operator([Outer(up, down'), Outer(down, up')], [1, 1])

# Or construct from arithmetic
σ_x = up * down' + down * up'
```

See also: [`Outer`](@ref), [`Identity`](@ref)
"""
struct Operator{S<:AbstractSpace} <: AbstractOperator{S}
    outers::Vector{Outer{S}}
    weights::Vector
    
    function Operator{S}(outers::Vector{Outer{S}}, weights::Vector) where S
        length(outers) == length(weights) || throw(ArgumentError("outers and weights must have same length"))
        isempty(outers) && throw(ArgumentError("Operator must have at least one outer product"))
        new{S}(outers, weights)
    end
    
    # Also allow calling without type parameter
    function Operator(outers::Vector{Outer{S}}, weights::Vector) where S
        length(outers) == length(weights) || throw(ArgumentError("outers and weights must have same length"))
        isempty(outers) && throw(ArgumentError("Operator must have at least one outer product"))
        new{S}(outers, weights)
    end
end

# Get space of an Operator
space(::Operator{S}) where S = S

# Display
function Base.show(io::IO, op::Operator)
    for (i, (outer, w)) in enumerate(zip(op.outers, op.weights))
        i > 1 && print(io, " + ")
        if !(w isa AbstractSymbolic) && isequal(w, 1)
            print(io, outer)
        else
            print(io, w, "·(", outer, ")")
        end
    end
end

# ============== Identity Operator ==============

@doc """
    Identity{S<:AbstractSpace}(space)

The identity operator on space S. Basis-independent.

# Examples
```julia
H, Hb = HilbertSpace(:H, 2)
I = Identity(H)

ψ = Ket(Hb, :ψ)
I * ψ  # → |ψ⟩
```

See also: [`Operator`](@ref), [`Outer`](@ref)
"""
struct Identity{S<:AbstractSpace} <: AbstractOperator{S}
    space::S
    
    function Identity(s::S) where S<:AbstractSpace
        new{S}(s)
    end
end

# Display
Base.show(io::IO, ::Identity) = print(io, "𝕀")

# ============== Function-based Operator ==============

@doc """
    FunctionOperator{S<:AbstractSpace, B<:AbstractBasis}(basis, action; adjoint_action=nothing, name=:F)
    FunctionOperator(basis) do ket ... end

Function-based operator that applies a user-defined function to kets in a specific basis.
The function receives a ket in the operator's basis and returns any AbstractKet.

If the input ket is in a different basis, the operator automatically applies a basis
transform before applying the action (if a transform is defined).

Constructed with do-block syntax:
```julia
F_op = FunctionOperator(basis) do ket
    # Process ket in 'basis', return AbstractKet or Number
    ...
end
```

Optionally provide `adjoint_action` for the adjoint operator and `name` for display.

# Examples
```julia
F, Fb = HilbertSpace(:Fock, Inf)

# Annihilation operator: â|n⟩ = √n |n-1⟩
â = FunctionOperator(Fb) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end

# With adjoint (creation operator):
â = FunctionOperator(Fb; 
    adjoint_action = ket -> begin
        n = parse(Int, string(ket.index))
        √(n + 1) * Ket(Fb, n + 1)
    end,
    name = :â
) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end

# Cross-basis: if ket is in different basis, transform is applied first
Fb2 = Basis(F, :energy)
# Define transform between bases first
define_transform!(typeof(Fb2), typeof(Fb)) do k
    # ... transformation logic ...
end
â * Ket(Fb2, :E0)  # Transforms to Fb basis first, then applies â
```

See also: [`Operator`](@ref), [`Outer`](@ref), [`define_transform!`](@ref)
"""
struct FunctionOperator{S<:AbstractSpace, B<:AbstractBasis} <: AbstractOperator{S}
    basis::B
    action::Function
    adjoint_action::Union{Function, Nothing}
    name::Symbol
    
    function FunctionOperator{S,B}(basis::B, action::Function, adjoint_action::Union{Function, Nothing}, name::Symbol) where {S<:AbstractSpace, B<:AbstractBasis}
        new{S,B}(basis, action, adjoint_action, name)
    end
end

# Constructor with do-block
function FunctionOperator(action::F, basis::B; adjoint_action::Union{Function, Nothing}=nothing, name::Symbol=:F) where {F<:Function, B<:AbstractBasis}
    S = space(basis)  # S is a type, not an instance
    FunctionOperator{S,B}(basis, action, adjoint_action, name)
end

# Regular constructor
function FunctionOperator(basis::B, action::F; adjoint_action::Union{Function, Nothing}=nothing, name::Symbol=:F) where {F<:Function, B<:AbstractBasis}
    S = space(basis)  # S is a type, not an instance
    FunctionOperator{S,B}(basis, action, adjoint_action, name)
end

Base.adjoint(op::FunctionOperator) = isnothing(op.adjoint_action) ? throw(ErrorException("Adjoint action not defined for this FunctionOperator")) : FunctionOperator(op.basis, op.adjoint_action, adjoint_action=op.action, name=Symbol(string(op.name), "'"))

# Display
Base.show(io::IO, op::FunctionOperator) = print(io, op.name)

# ============== OperatorSum: Lazy Sum of Operators ==============

@doc """
    OperatorSum{S<:AbstractSpace}

Lazy container for sum of operators. Preserves different operator types
(Identity, Operator, FunctionOperator, etc.) without forcing evaluation.

Enables expressions like:
- Op + c·𝕀 (operator plus scaled identity)
- Op₁ + Op₂ + FuncOp (mixing operator types)

Application evaluates lazily: (A + B)|ψ⟩ = A|ψ⟩ + B|ψ⟩

# Examples
```julia
H, Hb = HilbertSpace(:H, 2)
σ_z = Ket(Hb, :↑) * Ket(Hb, :↑)' - Ket(Hb, :↓) * Ket(Hb, :↓)'
H_shifted = σ_z + 2 * Identity(H)  # σ_z + 2𝕀

# Application
ψ = Ket(Hb, :↑)
H_shifted * ψ  # → σ_z|↑⟩ + 2|↑⟩ = 3|↑⟩
```

See also: [`Operator`](@ref), [`Identity`](@ref), [`Outer`](@ref)
"""
struct OperatorSum{S<:AbstractSpace} <: AbstractOperator{S}
    operators::Vector{<:AbstractOperator{S}}
    weights::Vector{<:Number}
    
    function OperatorSum(ops::Vector{<:AbstractOperator{S}}, weights::Vector{<:Number}) where S
        length(ops) == length(weights) || throw(ArgumentError("operators and weights must have same length"))
        isempty(ops) && throw(ArgumentError("OperatorSum requires at least one operator"))
        
        # All operators must act on same space (S is already constrained by type)
        new{S}(ops, weights)
    end
end

# Constructor: single operator + weight
OperatorSum(op::AbstractOperator{S}, weight::Number) where S = OperatorSum([op], [weight])

# Show method
function Base.show(io::IO, opsum::OperatorSum{S}) where S
    print(io, "OperatorSum on $(S): ")
    for (i, (op, w)) in enumerate(zip(opsum.operators, opsum.weights))
        i > 1 && print(io, " + ")
        !isone(w) && print(io, "($w)×")
        print(io, typeof(op).name.name)
    end
end

# ============== Trace ==============

@doc """
    tr(op::AbstractOperator)

Compute the trace of an operator.

For an outer product |ψ⟩⟨ϕ|, the trace is ⟨ϕ|ψ⟩.
For a sum of weighted outers Σᵢ wᵢ|ψᵢ⟩⟨ϕᵢ|, the trace is Σᵢ wᵢ⟨ϕᵢ|ψᵢ⟩.

# Examples
```julia
H, Hb = HilbertSpace(:H, 2)
up = Ket(Hb, :↑)
down = Ket(Hb, :↓)

# Projector - trace is 1
P_up = up * up'
tr(P_up)  # → 1

# Off-diagonal - trace is 0 (orthonormal basis)
σ_plus = up * down'
tr(σ_plus)  # → 0

# σ_x = |↑⟩⟨↓| + |↓⟩⟨↑| - trace is 0
σ_x = up * down' + down * up'
tr(σ_x)  # → 0
```

See also: [`Outer`](@ref), [`Operator`](@ref)
"""
function tr end

# Trace of Outer: ⟨ϕ|ψ⟩
tr(op::Outer) = op.bra * op.ket

# Trace of Operator: collect all kets and bras, sum over all pairs
# tr(Σᵢ wᵢ|ψᵢ⟩⟨ϕᵢ|) = Σᵢⱼ wᵢ ⟨ϕⱼ|ψᵢ⟩
function tr(op::Operator)
    # Collect all weighted kets and all bras from the outers
    weighted_kets = [(outer.ket, w) for (outer, w) in zip(op.outers, op.weights)]
    bras = [outer.bra for outer in op.outers]
    
    # Sum over all (bra, ket) pairs
    result = 0
    for bra in bras
        for (ket, w) in weighted_kets
            result = simplify(result + w * (bra * ket))
        end
    end
    return result
end

# Trace of Identity: dimension of the space
function tr(op::Identity{S}) where S<:AbstractSpace{name, dims} where {name, dims}
    # dims is a tuple of dimensions (e.g., (2,) or (2, 3) or (nothing,))
    # Total dimension is the product, or infinite if any is nothing
    if any(isnothing, dims)
        return Sym(:∞)  # Infinite-dimensional space
    end
    return prod(dims)
end

# Trace of OperatorSum: sum of traces
function tr(opsum::OperatorSum)
    result = 0
    for (op, w) in zip(opsum.operators, opsum.weights)
        result = simplify(result + w * tr(op))
    end
    return result
end
