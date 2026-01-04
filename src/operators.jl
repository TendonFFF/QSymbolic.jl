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

# ============== Operator-Operator Multiplication ==============

@doc """
Operator multiplication using bra-ket contraction.

For A = Σᵢ wᵢ|aᵢ⟩⟨bᵢ| and B = Σⱼ vⱼ|cⱼ⟩⟨dⱼ|:
    A * B = Σᵢⱼ wᵢvⱼ⟨bᵢ|cⱼ⟩ |aᵢ⟩⟨dⱼ|

The contraction ⟨bᵢ|cⱼ⟩ uses bra-ket arithmetic, which handles:
- Same basis, same index: 1
- Same basis, different index: 0
- Cross-basis: InnerProduct or custom transform
"""
function Base.:*(op1::Outer{S}, op2::Outer{S}) where S
    # (|a⟩⟨b|) * (|c⟩⟨d|) = ⟨b|c⟩ |a⟩⟨d|
    contraction = op1.bra * op2.ket
    
    if iszero(contraction)
        # Result is zero operator - return a zero-weighted outer product
        return 0 * Outer(op1.ket, op2.bra)
    else
        # Result is contraction * |a⟩⟨d|
        return contraction * Outer(op1.ket, op2.bra)
    end
end

function Base.:*(op1::Operator{S}, op2::Outer{S}) where S
    # Σᵢ wᵢ|aᵢ⟩⟨bᵢ| * |c⟩⟨d| = Σᵢ wᵢ⟨bᵢ|c⟩ |aᵢ⟩⟨d|
    result_outers = Outer{S}[]
    result_weights = []
    
    for (outer1, w1) in zip(op1.outers, op1.weights)
        contraction = outer1.bra * op2.ket
        if !iszero(contraction)
            push!(result_outers, Outer(outer1.ket, op2.bra))
            push!(result_weights, w1 * contraction)
        end
    end
    
    isempty(result_outers) && return 0
    return Operator{S}(result_outers, result_weights)
end

function Base.:*(op1::Outer{S}, op2::Operator{S}) where S
    # |a⟩⟨b| * Σⱼ vⱼ|cⱼ⟩⟨dⱼ| = Σⱼ vⱼ⟨b|cⱼ⟩ |a⟩⟨dⱼ|
    result_outers = Outer{S}[]
    result_weights = []
    
    for (outer2, w2) in zip(op2.outers, op2.weights)
        contraction = op1.bra * outer2.ket
        if !iszero(contraction)
            push!(result_outers, Outer(op1.ket, outer2.bra))
            push!(result_weights, w2 * contraction)
        end
    end
    
    isempty(result_outers) && return 0
    return Operator{S}(result_outers, result_weights)
end

function Base.:*(op1::Operator{S}, op2::Operator{S}) where S
    # Σᵢ wᵢ|aᵢ⟩⟨bᵢ| * Σⱼ vⱼ|cⱼ⟩⟨dⱼ| = Σᵢⱼ wᵢvⱼ⟨bᵢ|cⱼ⟩ |aᵢ⟩⟨dⱼ|
    result_outers = Outer{S}[]
    result_weights = []
    
    for (outer1, w1) in zip(op1.outers, op1.weights)
        for (outer2, w2) in zip(op2.outers, op2.weights)
            contraction = outer1.bra * outer2.ket
            if !iszero(contraction)
                push!(result_outers, Outer(outer1.ket, outer2.bra))
                push!(result_weights, w1 * w2 * contraction)
            end
        end
    end
    
    isempty(result_outers) && return 0
    return Operator{S}(result_outers, result_weights)
end

# Identity operator multiplication
Base.:*(::Identity{S}, op::AbstractOperator{S}) where S = op
Base.:*(op::AbstractOperator{S}, ::Identity{S}) where S = op

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
""" FunctionOperator
struct FunctionOperator{S<:AbstractSpace, B<:AbstractBasis} <: AbstractOperator{S}
    space::S
    basis::B
    action::Function
    adjoint_action::Union{Function, Nothing}
    name::Symbol
    
    function FunctionOperator{S,B}(space::S, basis::B, action::Function, adjoint_action::Union{Function, Nothing}, name::Symbol) where {S<:AbstractSpace, B<:AbstractBasis}
        space(basis) == space || throw(ArgumentError("Basis must be in the same space as operator"))
        new{S,B}(space, basis, action, adjoint_action, name)
    end
end

# Constructor with do-block
function FunctionOperator(action::F, basis::B; adjoint_action::Union{Function, Nothing}=nothing, name::Symbol=:F) where {F<:Function, B<:AbstractBasis}
    s = space(basis)
    FunctionOperator{typeof(s),typeof(basis)}(s, basis, action, adjoint_action, name)
end

# Regular constructor
function FunctionOperator(basis::B, action::F; adjoint_action::Union{Function, Nothing}=nothing, name::Symbol=:F) where {F<:Function, B<:AbstractBasis}
    s = space(basis)
    FunctionOperator{typeof(s),typeof(basis)}(s, basis, action, adjoint_action, name)
end

# Apply function operator with automatic basis transform
function Base.:*(op::FunctionOperator, ket::AbstractKet)
    ket_basis = basis(ket)
    
    # If ket is in different basis, try to transform
    if ket_basis != op.basis
        # Check if transformation is defined
        if has_transform(typeof(ket_basis), typeof(op.basis))
            # Apply transform to operator's basis
            ket = transform(ket, op.basis)
        else
            throw(ArgumentError("No transform defined from basis $(ket_basis) to $(op.basis). Define one with define_transform!"))
        end
    end
    
    # Now apply the action (ket is in correct basis)
    op.action(ket)
end

# Make callable
(op::FunctionOperator)(ket::AbstractKet) = op * ket

# Display
Base.show(io::IO, op::FunctionOperator) = print(io, op.name)

# Adjoint (wrapper)
struct AdjointFunctionOperator{S<:AbstractSpace, B<:AbstractBasis} <: AbstractOperator{S}
    parent::FunctionOperator{S,B}
end

Base.adjoint(op::FunctionOperator) = AdjointFunctionOperator(op)
Base.adjoint(op::AdjointFunctionOperator) = op.parent

Base.show(io::IO, op::AdjointFunctionOperator) = print(io, op.parent.name, "†")

# Apply adjoint with automatic basis transform
function Base.:*(op::AdjointFunctionOperator, ket::AbstractKet)
    if !isnothing(op.parent.adjoint_action)
        ket_basis = basis(ket)
        
        # If ket is in different basis, try to transform
        if ket_basis != op.parent.basis
            # Check if transformation is defined
            if has_transform(typeof(ket_basis), typeof(op.parent.basis))
                # Apply transform to operator's basis
                ket = transform(ket, op.parent.basis)
            else
                throw(ArgumentError("No transform defined from basis $(ket_basis) to $(op.parent.basis). Define one with define_transform!"))
            end
        end
        
        # Apply adjoint action
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
""" OperatorSum
struct OperatorSum{S<:AbstractSpace} <: AbstractOperator{S}
    operators::Vector{AbstractOperator{S}}
    weights::Vector{Number}
    space::S
    
    function OperatorSum(ops::Vector{<:AbstractOperator{S}}, weights::Vector{<:Number}) where S
        length(ops) == length(weights) || throw(ArgumentError("operators and weights must have same length"))
        isempty(ops) && throw(ArgumentError("OperatorSum requires at least one operator"))
        
        # All operators must act on same space
        sp = space(ops[1])
        all(space(op) == sp for op in ops) || throw(ArgumentError("All operators must act on same space"))
        
        new{S}(ops, weights, sp)
    end
end

# Constructor: single operator + weight
OperatorSum(op::AbstractOperator{S}, weight::Number) where S = OperatorSum([op], [weight])

# Show method
function Base.show(io::IO, opsum::OperatorSum)
    print(io, "OperatorSum on $(opsum.space): ")
    for (i, (op, w)) in enumerate(zip(opsum.operators, opsum.weights))
        i > 1 && print(io, " + ")
        !isone(w) && print(io, "($w)×")
        print(io, typeof(op).name.name)
    end
end

# ============== OperatorSum Arithmetic ==============

# Addition: OperatorSum + anything
Base.:+(opsum::OperatorSum{S}, op::AbstractOperator{S}) where S = 
    OperatorSum(vcat(opsum.operators, [op]), vcat(opsum.weights, [1]))

Base.:+(opsum::OperatorSum{S}, op::WeightedOperator{S}) where S = 
    OperatorSum(vcat(opsum.operators, [op.operator]), vcat(opsum.weights, [op.weight]))

Base.:+(op::AbstractOperator{S}, opsum::OperatorSum{S}) where S = opsum + op

# Addition: OperatorSum + OperatorSum
Base.:+(opsum1::OperatorSum{S}, opsum2::OperatorSum{S}) where S =
    OperatorSum(vcat(opsum1.operators, opsum2.operators), vcat(opsum1.weights, opsum2.weights))

# Addition: Operator + Identity (the key use case!)
Base.:+(op::Operator{S}, id::Identity{S}) where S = OperatorSum([op, id], [1, 1])
Base.:+(id::Identity{S}, op::Operator{S}) where S = OperatorSum([id, op], [1, 1])

Base.:+(op::Outer{S}, id::Identity{S}) where S = OperatorSum([op, id], [1, 1])
Base.:+(id::Identity{S}, op::Outer{S}) where S = OperatorSum([id, op], [1, 1])

Base.:+(op::FunctionOperator{S}, id::Identity{S}) where S = OperatorSum([op, id], [1, 1])
Base.:+(id::Identity{S}, op::FunctionOperator{S}) where S = OperatorSum([id, op], [1, 1])

# Scalar multiplication
Base.:*(c::Number, opsum::OperatorSum) = OperatorSum(opsum.operators, c .* opsum.weights)
Base.:*(opsum::OperatorSum, c::Number) = c * opsum

# ============== OperatorSum Application ==============

# OperatorSum acting on ket: (Σᵢ wᵢ Opᵢ)|ψ⟩ = Σᵢ wᵢ(Opᵢ|ψ⟩)
function Base.:*(opsum::OperatorSum{S}, ket::AbstractKet) where S
    total = nothing
    for (op, w) in zip(opsum.operators, opsum.weights)
        result = op * ket
        if !iszero(result)
            term = w * result
            total = isnothing(total) ? term : total + term
        end
    end
    isnothing(total) ? 0 : total
end

# Helper type for weighted operators (used internally)
struct WeightedOperator{S<:AbstractSpace} <: AbstractOperator{S}
    operator::AbstractOperator{S}
    weight::Number
end
Base.:*(w::Number, op::AbstractOperator) = WeightedOperator(op, w)
Base.:*(op::AbstractOperator, w::Number) = w * op
Base.:*(wop::WeightedOperator, ket::AbstractKet) = wop.weight * (wop.operator * ket)

# ============== Bra-Operator Contractions ==============

@doc """
Bra-operator contractions: ⟨ψ|Â

For most operators, this is computed as (Â†|ψ⟩)† = ⟨ψ|Â
This uses the adjoint of the operator to act on the ket, then takes the adjoint of the result.
"""

# Bra * Outer: ⟨χ|(|ψ⟩⟨ϕ|) = ⟨χ|ψ⟩⟨ϕ| 
function Base.:*(bra::AbstractBra, op::Outer)
    # ⟨χ|(|ψ⟩⟨ϕ|) = (|ϕ⟩⟨ψ|)|χ⟩)† = (⟨ψ|χ⟩|ϕ⟩)† = ⟨ψ|χ⟩⟨ϕ|
    inner = adjoint(op.ket) * bra
    if iszero(inner)
        return 0
    else
        return inner * op.bra
    end
end

# Bra * Operator: ⟨χ|(Σᵢ wᵢ|ψᵢ⟩⟨ϕᵢ|) = Σᵢ wᵢ⟨χ|ψᵢ⟩⟨ϕᵢ|
function Base.:*(bra::AbstractBra, op::Operator)
    result = nothing
    for (outer, w) in zip(op.outers, op.weights)
        term = bra * outer
        if !iszero(term)
            weighted_term = w * term
            result = isnothing(result) ? weighted_term : result + weighted_term
        end
    end
    isnothing(result) ? 0 : result
end

# Bra * Identity: ⟨ψ|𝕀 = ⟨ψ|
Base.:*(bra::AbstractBra, ::Identity) = bra

# Bra * FunctionOperator: ⟨ψ|F̂
function Base.:*(bra::AbstractBra, op::FunctionOperator)
    # Try to use adjoint action if available
    if !isnothing(op.adjoint_action)
        # Apply adjoint operator: (F̂†|ψ⟩)† 
        ket = adjoint(bra)
        result_ket = op' * ket
        return adjoint(result_ket)
    else
        throw(ArgumentError("Bra-operator contraction requires adjoint_action to be defined for FunctionOperator $(op.name)"))
    end
end

# Bra * AdjointFunctionOperator: ⟨ψ|F̂†
function Base.:*(bra::AbstractBra, op::AdjointFunctionOperator)
    # Apply the original operator's action to the adjoint of the bra
    ket = adjoint(bra)
    result_ket = op.parent * ket
    return adjoint(result_ket)
end

# Bra * OperatorSum: ⟨ψ|(Σᵢ wᵢ Opᵢ) = Σᵢ wᵢ⟨ψ|Opᵢ
function Base.:*(bra::AbstractBra, opsum::OperatorSum)
    result = nothing
    for (op, w) in zip(opsum.operators, opsum.weights)
        term = bra * op
        if !iszero(term)
            weighted_term = w * term
            result = isnothing(result) ? weighted_term : result + weighted_term
        end
    end
    isnothing(result) ? 0 : result
end

# Bra * WeightedOperator
Base.:*(bra::AbstractBra, wop::WeightedOperator) = wop.weight * (bra * wop.operator)
