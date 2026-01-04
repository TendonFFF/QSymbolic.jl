# Quantum operators - Restructured to work with bra-ket arithmetic

# ============== Abstract Operator Type ==============

@doc """
    AbstractOperator{S<:AbstractSpace}

Abstract supertype for all quantum operators acting on space `S`.
Operators use bra-ket arithmetic for contraction.
""" AbstractOperator
abstract type AbstractOperator{S<:AbstractSpace} end

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
""" Outer
struct Outer{S<:AbstractSpace} <: AbstractOperator{S}
    ket::AbstractKet
    bra::AbstractBra
    space::S
    
    function Outer(ket::AbstractKet, bra::AbstractBra)
        # Verify same space
        ket_space = space(basis(ket))
        bra_space = space(basis(bra))
        ket_space == bra_space || throw(ArgumentError("Ket and bra must be in the same space"))
        new{typeof(ket_space)}(ket, bra, ket_space)
    end
end

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

# Apply outer product to ket: (|ψ⟩⟨ϕ|)|χ⟩ = ⟨ϕ|χ⟩ |ψ⟩
function Base.:*(op::Outer, ket::AbstractKet)
    # Use bra-ket contraction
    inner = op.bra * ket
    if !(inner isa AbstractSymbolic) && iszero(inner)
        return 0
    else
        return inner * op.ket
    end
end

# Adjoint: (|ψ⟩⟨ϕ|)† = |ϕ⟩⟨ψ|
Base.adjoint(op::Outer) = Outer(adjoint(op.bra), adjoint(op.ket))

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
""" Operator
struct Operator{S<:AbstractSpace} <: AbstractOperator{S}
    outers::Vector{Outer{S}}
    weights::Vector
    space::S
    
    function Operator(outers::Vector{Outer{S}}, weights::Vector) where S
        length(outers) == length(weights) || throw(ArgumentError("outers and weights must have same length"))
        isempty(outers) && throw(ArgumentError("Operator must have at least one outer product"))
        space_val = space(outers[1])
        new{S}(outers, weights, space_val)
    end
end

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

# Apply operator to ket: decomposes to bra-ket contractions
function Base.:*(op::Operator, ket::AbstractKet)
    total = nothing
    for (outer, w) in zip(op.outers, op.weights)
        result = outer * ket
        if !iszero(result)
            term = w * result
            total = isnothing(total) ? term : total + term
        end
    end
    isnothing(total) ? 0 : total
end

# Adjoint: (Σᵢ wᵢ |ψᵢ⟩⟨ϕᵢ|)† = Σᵢ wᵢ* |ϕᵢ⟩⟨ψᵢ|
function Base.adjoint(op::Operator{S}) where S
    Operator{S}([adjoint(o) for o in op.outers], [conj(w) for w in op.weights])
end

# Addition of operators
function Base.:+(op1::Operator{S}, op2::Operator{S}) where S
    Operator{S}(vcat(op1.outers, op2.outers), vcat(op1.weights, op2.weights))
end

function Base.:+(outer1::Outer{S}, outer2::Outer{S}) where S
    Operator{S}([outer1, outer2], [1, 1])
end

function Base.:+(op::Operator{S}, outer::Outer{S}) where S
    Operator{S}(vcat(op.outers, [outer]), vcat(op.weights, [1]))
end

function Base.:+(outer::Outer{S}, op::Operator{S}) where S
    Operator{S}(vcat([outer], op.outers), vcat([1], op.weights))
end

# Subtraction
Base.:-(op1::AbstractOperator{S}, op2::AbstractOperator{S}) where S = op1 + (-1 * op2)

# Scalar multiplication
function Base.:*(a::Number, op::Operator{S}) where S
    Operator{S}(op.outers, [a * w for w in op.weights])
end

function Base.:*(a::Number, outer::Outer{S}) where S
    Operator{S}([outer], [a])
end

Base.:*(op::AbstractOperator, a::Number) = a * op
Base.:/(op::AbstractOperator, a::Number) = (1/a) * op

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
""" Identity
struct Identity{S<:AbstractSpace} <: AbstractOperator{S}
    space::S
    
    function Identity(s::S) where S<:AbstractSpace
        new{S}(s)
    end
end

# Apply identity
Base.:*(::Identity, ket::AbstractKet) = ket

# Operator algebra with identity
Base.:*(::Identity{S}, op::AbstractOperator{S}) where S = op
Base.:*(op::AbstractOperator{S}, ::Identity{S}) where S = op
Base.:*(::Identity{S}, ::Identity{S}) where S = Identity(S)

# Adjoint
Base.adjoint(id::Identity) = id

# Display
Base.show(io::IO, ::Identity) = print(io, "𝕀")

# ============== Function-based Operator ==============

@doc """
    FunctionOperator{S<:AbstractSpace}(space, action; adjoint_action=nothing)

An operator defined by a function that maps kets to kets (or numbers).
The function can return any AbstractKet, regardless of basis.

Constructed with do-block syntax:
```julia
F_op = FunctionOperator(space) do ket
    # Process ket, return AbstractKet or Number
    ...
end
```

Optionally provide `adjoint_action` for the adjoint operator.

# Examples
```julia
F, Fb = FockSpace(:F)

# Annihilation operator: â|n⟩ = √n |n-1⟩
â = FunctionOperator(F) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end

# With adjoint (creation operator):
â = FunctionOperator(F; 
    adjoint_action = ket -> begin
        n = parse(Int, string(ket.index))
        √(n + 1) * Ket(Fb, n + 1)
    end
) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end
```

See also: [`Operator`](@ref), [`Outer`](@ref)
""" FunctionOperator
struct FunctionOperator{S<:AbstractSpace} <: AbstractOperator{S}
    space::S
    action::Function
    adjoint_action::Union{Function, Nothing}
    name::Symbol
    
    function FunctionOperator{S}(space::S, action::Function, adjoint_action::Union{Function, Nothing}, name::Symbol) where S<:AbstractSpace
        new{S}(space, action, adjoint_action, name)
    end
end

# Constructor with do-block
function FunctionOperator(action::F, space::S; adjoint_action::Union{Function, Nothing}=nothing, name::Symbol=:F) where {F<:Function, S<:AbstractSpace}
    FunctionOperator{S}(space, action, adjoint_action, name)
end

# Regular constructor
function FunctionOperator(space::S, action::F; adjoint_action::Union{Function, Nothing}=nothing, name::Symbol=:F) where {F<:Function, S<:AbstractSpace}
    FunctionOperator{S}(space, action, adjoint_action, name)
end

# Apply function operator
Base.:*(op::FunctionOperator, ket::AbstractKet) = op.action(ket)

# Make callable
(op::FunctionOperator)(ket::AbstractKet) = op * ket

# Display
Base.show(io::IO, op::FunctionOperator) = print(io, op.name)

# Adjoint (wrapper)
struct AdjointFunctionOperator{S<:AbstractSpace} <: AbstractOperator{S}
    parent::FunctionOperator{S}
end

Base.adjoint(op::FunctionOperator) = AdjointFunctionOperator(op)
Base.adjoint(op::AdjointFunctionOperator) = op.parent

Base.show(io::IO, op::AdjointFunctionOperator) = print(io, op.parent.name, "†")

# Apply adjoint
function Base.:*(op::AdjointFunctionOperator, ket::AbstractKet)
    if !isnothing(op.parent.adjoint_action)
        op.parent.adjoint_action(ket)
    else
        throw(ArgumentError("Adjoint action not defined for $(op.parent.name)"))
    end
end

# ============== Partial Application on ProductKet ==============

@doc """
Outer{S1} acting on ProductKet with space S1⊗S2:
Applies the operator to the S1 component, leaves S2 unchanged.

Result: |result_S1⟩ ⊗ |original_S2⟩

# Example
```julia
H1, Hb1 = HilbertSpace(:A, 2)
H2, Hb2 = HilbertSpace(:B, 2)

up1 = Ket(Hb1, :↑)
down1 = Ket(Hb1, :↓)
ψ2 = Ket(Hb2, :ψ)

# Operator on H1
P_up = up1 * up1'

# Product state
state = up1 ⊗ ψ2

# P_up acts only on H1 component
P_up * state  # → |↑⟩ ⊗ |ψ⟩
```
"""
function Base.:*(op::Outer{S1}, ket::ProductKet) where S1
    # Check if operator's space matches any component of the product ket
    # Find which component(s) match
    matched_indices = Int[]
    for (i, k) in enumerate(ket.kets)
        if space(basis(k)) == op.space
            push!(matched_indices, i)
        end
    end
    
    if isempty(matched_indices)
        throw(ArgumentError("Operator space $(op.space) doesn't match any component of ProductKet"))
    end
    
    if length(matched_indices) > 1
        throw(ArgumentError("Operator space matches multiple components. Specify which one."))
    end
    
    idx = matched_indices[1]
    
    # Apply operator to the matched component
    result_component = op * ket.kets[idx]
    
    # Handle zero result
    iszero(result_component) && return 0
    
    # Reconstruct ProductKet with new component
    new_kets = copy(ket.kets)
    
    if result_component isa Ket
        new_kets[idx] = result_component
        return ProductKet(new_kets)
    elseif result_component isa WeightedKet
        new_kets[idx] = result_component.ket
        return result_component.weight * ProductKet(new_kets)
    elseif result_component isa SumKet
        # Distribute: (Σᵢ wᵢ|ψᵢ⟩) ⊗ |other⟩ = Σᵢ wᵢ(|ψᵢ⟩ ⊗ |other⟩)
        result_kets = ProductKet[]
        result_weights = []
        for (k, w) in zip(result_component.kets, result_component.weights)
            temp_kets = copy(new_kets)
            temp_kets[idx] = k
            push!(result_kets, ProductKet(temp_kets))
            push!(result_weights, w)
        end
        return SumKet(result_kets, result_weights)
    else
        # Number or other type
        temp_kets = copy(new_kets)
        return result_component * ProductKet(temp_kets)
    end
end

# Operator container acting on ProductKet
function Base.:*(op::Operator{S1}, ket::ProductKet) where S1
    total = nothing
    for (outer, w) in zip(op.outers, op.weights)
        result = outer * ket
        if !iszero(result)
            term = w * result
            total = isnothing(total) ? term : total + term
        end
    end
    isnothing(total) ? 0 : total
end

# FunctionOperator acting on ProductKet
function Base.:*(op::FunctionOperator{S1}, ket::ProductKet) where S1
    # Check if operator's space matches product ket's composite space
    composite_space = space(basis(ket))
    if composite_space == op.space
        # Operator acts on entire product ket
        return op.action(ket)
    end
    
    # Otherwise, try partial application (if space matches a component)
    matched_indices = Int[]
    for (i, k) in enumerate(ket.kets)
        if space(basis(k)) == op.space
            push!(matched_indices, i)
        end
    end
    
    if length(matched_indices) == 1
        # Apply to single component
        idx = matched_indices[1]
        result_component = op * ket.kets[idx]
        iszero(result_component) && return 0
        
        new_kets = copy(ket.kets)
        if result_component isa Ket
            new_kets[idx] = result_component
            return ProductKet(new_kets)
        elseif result_component isa WeightedKet
            new_kets[idx] = result_component.ket
            return result_component.weight * ProductKet(new_kets)
        elseif result_component isa SumKet
            result_kets = ProductKet[]
            result_weights = []
            for (k, w) in zip(result_component.kets, result_component.weights)
                temp_kets = copy(new_kets)
                temp_kets[idx] = k
                push!(result_kets, ProductKet(temp_kets))
                push!(result_weights, w)
            end
            return SumKet(result_kets, result_weights)
        else
            return result_component * ProductKet(new_kets)
        end
    else
        throw(ArgumentError("Cannot apply FunctionOperator: space mismatch"))
    end
end

# ============== Zero checking ==============

Base.iszero(x::Number) = !(x isa AbstractSymbolic) && x == 0
Base.iszero(::AbstractKet) = false
Base.iszero(::AbstractOperator) = false

# ============== Tensor Product of Operators ==============

@doc """
    op1 ⊗ op2

Tensor product of two operators: Â ⊗ B̂. 
Acts on composite states as: (Â ⊗ B̂)(|ψ⟩ ⊗ |ϕ⟩) = (Â|ψ⟩) ⊗ (B̂|ϕ⟩)
"""
function ⊗(op1::AbstractOperator{S1}, op2::AbstractOperator{S2}) where {S1, S2}
    # Create a combined operator that acts on S1⊗S2
    combined_space = S1 ⊗ S2  # Assumes spaces have ⊗ method
    
    FunctionOperator(combined_space) do ket
        if ket isa ProductKet && length(ket.kets) >= 2
            # Find components matching S1 and S2
            # This is simplified - full implementation would need better matching
            result1 = op1 * ket.kets[1]
            result2 = op2 * ket.kets[2]
            
            (iszero(result1) || iszero(result2)) && return 0
            
            # Combine results
            if result1 isa Ket && result2 isa Ket
                return ProductKet([result1, result2])
            else
                # Handle weighted/sum combinations
                return result1 ⊗ result2
            end
        else
            throw(ArgumentError("Tensor product operators require ProductKet"))
        end
    end
end

# ============== Convenience constructors ==============

@doc """
    create_ladder_operators(space, basis; name=:a)

Create bosonic annihilation and creation operators for a Fock space.
Returns `(â, â†)` where:
- â|n⟩ = √n |n-1⟩  (annihilation)
- â†|n⟩ = √(n+1) |n+1⟩  (creation)
"""
function create_ladder_operators(space::S, basis::B; name::Symbol=:a) where {S<:AbstractSpace, B<:AbstractBasis}
    â = FunctionOperator(space; name=name,
        adjoint_action = ket -> begin
            n = ket.index isa Symbol ? parse(Int, string(ket.index)) : Int(ket.index)
            √(n + 1) * Ket(basis, n + 1)
        end
    ) do ket
        n = ket.index isa Symbol ? parse(Int, string(ket.index)) : Int(ket.index)
        n == 0 ? 0 : √n * Ket(basis, n - 1)
    end
    
    return â, â'
end
