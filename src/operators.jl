# Quantum operators

# ============== Abstract Operator Type ==============

@doc """
    AbstractOperator{B}

Abstract supertype for all quantum operators in basis `B`.
""" AbstractOperator
abstract type AbstractOperator{B<:AbstractBasis} end

# ============== Outer Product Operator |ψ⟩⟨ϕ| ==============

@doc """
    Operator{B}(ket, bra)
    ket * bra'  (automatically creates Operator)

An operator formed as the outer product |ψ⟩⟨ϕ|. This is the standard
quantum mechanical representation of operators.

When applied to a state: (|ψ⟩⟨ϕ|)|χ⟩ = |ψ⟩⟨ϕ|χ⟩ = ⟨ϕ|χ⟩ |ψ⟩

# Examples
```julia
H = HilbertSpace(:H, 2)
Zb = Basis(H, :z)
up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

# Projector onto |↑⟩
P_up = up * up'  # |↑⟩⟨↑|

# Apply: P_up|↑⟩ = |↑⟩, P_up|↓⟩ = 0
P_up * up   # → |↑⟩
P_up * down # → 0

# Ladder operator |↑⟩⟨↓|
σ_plus = up * down'
σ_plus * down  # → |↑⟩
```

See also: [`SumOperator`](@ref), [`ScaledOperator`](@ref), [`FunctionOperator`](@ref)
""" Operator
struct Operator{B<:AbstractBasis, T} <: AbstractOperator{B}
    ket::AbstractKet{B}
    bra::AbstractBra{B}
    coeff::T
    
    function Operator(ket::AbstractKet{B}, bra::AbstractBra{B}, coeff::T=1) where {B<:AbstractBasis, T}
        new{B, T}(ket, bra, coeff)
    end
end

# basis() - documented in basis.jl (additional method for operators)
basis(::AbstractOperator{B}) where B = B

# space() - documented in basis.jl (additional method for operators)
space(::AbstractOperator{B}) where B = space(B)

# Create operator from ket * bra
function Base.:*(ket::Ket{B}, bra::Bra{B}) where B<:AbstractBasis
    Operator(ket, bra, 1)
end

function Base.:*(ket::WeightedKet{B}, bra::Bra{B}) where B
    Operator(ket.ket, bra, ket.weight)
end

function Base.:*(ket::Ket{B}, bra::WeightedBra{B}) where B
    Operator(ket, bra.bra, bra.weight)
end

function Base.:*(ket::WeightedKet{B}, bra::WeightedBra{B}) where B
    Operator(ket.ket, bra.bra, ket.weight * bra.weight)
end

# Display
function Base.show(io::IO, op::Operator)
    if !(op.coeff isa AbstractSymbolic) && isequal(op.coeff, 1)
        print(io, _ket_str(op.ket), _bra_str(op.bra))
    else
        print(io, op.coeff, "·", _ket_str(op.ket), _bra_str(op.bra))
    end
end

# Helper functions for display
_ket_str(k::Ket) = "|$(k.index)⟩"
_bra_str(b::Bra) = "⟨$(b.index)|"

# Apply operator to ket: (|ψ⟩⟨ϕ|)|χ⟩ = ⟨ϕ|χ⟩ |ψ⟩
function Base.:*(op::Operator{B}, ket::Ket{B}) where B
    inner = Bra{B}(op.bra.index) * ket
    if !(inner isa AbstractSymbolic) && isequal(inner, 0)
        return 0
    else
        op.coeff * inner * op.ket
    end
end

function Base.:*(op::Operator{B}, ket::WeightedKet{B}) where B
    result = op * ket.ket
    result isa Number ? result * ket.weight : ket.weight * result
end

function Base.:*(op::Operator{B}, ket::SumKet{B,T}) where {B,T}
    total = nothing
    for (k, w) in zip(ket.kets, ket.weights)
        result = op * k
        term = result isa Number ? result * w : w * result
        if !iszero(term)
            total = isnothing(total) ? term : total + term
        end
    end
    isnothing(total) ? 0 : total
end

# Adjoint: (|ψ⟩⟨ϕ|)† = |ϕ⟩⟨ψ|
function Base.adjoint(op::Operator{B,T}) where {B,T}
    # Swap ket↔bra and conjugate coefficient
    new_ket = Ket{B}(op.bra.index)
    new_bra = Bra{B}(op.ket isa Ket ? op.ket.index : nothing)
    Operator(new_ket, new_bra, conj(op.coeff))
end

# ============== Sum of Operators ==============

@doc """
    SumOperator{B}

Sum of operators: Â + B̂. Created via operator addition.
""" SumOperator
struct SumOperator{B<:AbstractBasis} <: AbstractOperator{B}
    terms::Vector{<:AbstractOperator{B}}
end

# basis inherited from AbstractOperator

function Base.:+(op1::AbstractOperator{B}, op2::AbstractOperator{B}) where B
    terms1 = op1 isa SumOperator ? op1.terms : [op1]
    terms2 = op2 isa SumOperator ? op2.terms : [op2]
    SumOperator{B}(vcat(terms1, terms2))
end

function Base.:-(op1::AbstractOperator{B}, op2::AbstractOperator{B}) where B
    op1 + (-1 * op2)
end

function Base.show(io::IO, op::SumOperator)
    for (i, term) in enumerate(op.terms)
        i > 1 && print(io, " + ")
        print(io, term)
    end
end

# Apply sum of operators - specific ket types to avoid ambiguity
function Base.:*(op::SumOperator{B}, ket::Ket{B}) where B
    total = nothing
    for term in op.terms
        result = term * ket
        if !iszero(result)
            total = isnothing(total) ? result : total + result
        end
    end
    isnothing(total) ? 0 : total
end

function Base.:*(op::SumOperator{B}, ket::WeightedKet{B}) where B
    total = nothing
    for term in op.terms
        result = term * ket
        if !iszero(result)
            total = isnothing(total) ? result : total + result
        end
    end
    isnothing(total) ? 0 : total
end

function Base.:*(op::SumOperator{B}, ket::SumKet{B,T}) where {B,T}
    total = nothing
    for term in op.terms
        result = term * ket
        if !iszero(result)
            total = isnothing(total) ? result : total + result
        end
    end
    isnothing(total) ? 0 : total
end

# SumOperator on ProductKet (CompositeBasis)
function Base.:*(op::SumOperator{CompositeBasis{B1,B2}}, ket::ProductKet{B1,B2}) where {B1,B2}
    total = nothing
    for term in op.terms
        result = term * ket
        if !iszero(result)
            total = isnothing(total) ? result : total + result
        end
    end
    isnothing(total) ? 0 : total
end

# SumOperator on SumKet (CompositeBasis)
function Base.:*(op::SumOperator{CompositeBasis{B1,B2}}, ket::SumKet{B1,B2,T}) where {B1,B2,T}
    total = nothing
    for term in op.terms
        result = term * ket
        if !iszero(result)
            total = isnothing(total) ? result : total + result
        end
    end
    isnothing(total) ? 0 : total
end

function Base.adjoint(op::SumOperator{B}) where B
    SumOperator{B}([adjoint(t) for t in op.terms])
end

# ============== Scaled Operator ==============

@doc """
    ScaledOperator{B,T}

Operator multiplied by a scalar. Created via `scalar * operator`.
""" ScaledOperator
struct ScaledOperator{B<:AbstractBasis, T} <: AbstractOperator{B}
    coeff::T
    op::AbstractOperator{B}
end

# basis inherited from AbstractOperator

function Base.:*(a::Number, op::AbstractOperator{B}) where B
    ScaledOperator{B, typeof(a)}(a, op)
end

function Base.:*(op::AbstractOperator, a::Number)
    a * op
end

function Base.:/(op::AbstractOperator, a::Number)
    (1/a) * op
end

function Base.show(io::IO, op::ScaledOperator)
    print(io, op.coeff, "·(", op.op, ")")
end

function Base.:*(op::ScaledOperator{B}, ket::Ket{B}) where B
    result = op.op * ket
    result isa Number ? op.coeff * result : op.coeff * result
end

function Base.:*(op::ScaledOperator{B}, ket::WeightedKet{B}) where B
    result = op.op * ket
    result isa Number ? op.coeff * result : op.coeff * result
end

function Base.:*(op::ScaledOperator{B}, ket::SumKet{B,T}) where {B,T}
    result = op.op * ket
    result isa Number ? op.coeff * result : op.coeff * result
end

# ScaledOperator on ProductKet (CompositeBasis)
function Base.:*(op::ScaledOperator{CompositeBasis{B1,B2}}, ket::ProductKet{B1,B2}) where {B1,B2}
    result = op.op * ket
    result isa Number ? op.coeff * result : op.coeff * result
end

# ScaledOperator on SumKet (CompositeBasis)
function Base.:*(op::ScaledOperator{CompositeBasis{B1,B2}}, ket::SumKet{B1,B2,T}) where {B1,B2,T}
    result = op.op * ket
    result isa Number ? op.coeff * result : op.coeff * result
end

function Base.adjoint(op::ScaledOperator{B,T}) where {B,T}
    ScaledOperator{B,T}(conj(op.coeff), adjoint(op.op))
end

# Avoid double scaling
function Base.:*(a::Number, op::ScaledOperator{B,T}) where {B,T}
    ScaledOperator{B, promote_type(typeof(a),T)}(a * op.coeff, op.op)
end

# ============== Operator Product ==============

@doc """
    OperatorProduct{B}

Product of operators: ÂB̂. Created via operator multiplication.
""" OperatorProduct
struct OperatorProduct{B<:AbstractBasis} <: AbstractOperator{B}
    ops::Vector{<:AbstractOperator{B}}
end

# basis inherited from AbstractOperator

function Base.:*(op1::AbstractOperator{B}, op2::AbstractOperator{B}) where B
    ops1 = op1 isa OperatorProduct ? op1.ops : [op1]
    ops2 = op2 isa OperatorProduct ? op2.ops : [op2]
    OperatorProduct{B}(vcat(ops1, ops2))
end

function Base.show(io::IO, prod::OperatorProduct)
    for op in prod.ops
        print(io, "(", op, ")")
    end
end

# Apply product: (ÂB̂)|ψ⟩ = Â(B̂|ψ⟩) - specific ket types to avoid ambiguity
function Base.:*(prod::OperatorProduct{B}, ket::Ket{B}) where B
    result = ket
    for op in reverse(prod.ops)
        result = op * result
        iszero(result) && return 0
    end
    result
end

function Base.:*(prod::OperatorProduct{B}, ket::WeightedKet{B}) where B
    result = ket
    for op in reverse(prod.ops)
        result = op * result
        iszero(result) && return 0
    end
    result
end

function Base.:*(prod::OperatorProduct{B}, ket::SumKet{B,T}) where {B,T}
    result = ket
    for op in reverse(prod.ops)
        result = op * result
        iszero(result) && return 0
    end
    result
end

# OperatorProduct on ProductKet (CompositeBasis)
function Base.:*(prod::OperatorProduct{CompositeBasis{B1,B2}}, ket::ProductKet{B1,B2}) where {B1,B2}
    result = ket
    for op in reverse(prod.ops)
        result = op * result
        iszero(result) && return 0
    end
    result
end

# OperatorProduct on SumKet (CompositeBasis)
function Base.:*(prod::OperatorProduct{CompositeBasis{B1,B2}}, ket::SumKet{B1,B2,T}) where {B1,B2,T}
    result = ket
    for op in reverse(prod.ops)
        result = op * result
        iszero(result) && return 0
    end
    result
end

function Base.adjoint(prod::OperatorProduct{B}) where B
    OperatorProduct{B}(reverse([adjoint(op) for op in prod.ops]))
end

# ============== Function-based Operator ==============

@doc """
    FunctionOperator{B}(name, basis, action; adjoint_action=nothing)

A quantum operator defined by a function. Created via `FunctionOperator(name, basis) do ket ... end`.

The `action` function receives a `Ket` and should return:
- A number (scalar)
- A ket (`Ket`, `WeightedKet`, or `SumKet`)
- Zero (for annihilation of vacuum, etc.)

Optionally provide `adjoint_action` to define how the adjoint operator acts.
If not provided, applying the adjoint returns a symbolic `OpKet`.

Use this for operators that are easier to define procedurally (e.g., ladder operators).

# Example
```julia
F = FockSpace(:F)
Fb = Basis(F, :n)

# Annihilation operator: â|n⟩ = √n |n-1⟩
â = FunctionOperator(:â, Fb) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end

# With explicit adjoint action (creation operator):
â = FunctionOperator(:â, Fb; 
    adjoint_action = ket -> begin
        n = parse(Int, string(ket.index))
        √(n + 1) * Ket(Fb, n + 1)
    end
) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end
```

See also: [`Operator`](@ref), [`create_ladder_operators`](@ref)
""" FunctionOperator
struct FunctionOperator{B<:AbstractBasis} <: AbstractOperator{B}
    name::Symbol
    basis::B
    action::Function
    adjoint_action::Union{Function, Nothing}
    
    function FunctionOperator{B}(name::Symbol, basis::B, action::Function, adjoint_action::Union{Function, Nothing}=nothing) where B<:AbstractBasis
        new{B}(name, basis, action, adjoint_action)
    end
end

# basis inherited from AbstractOperator
# Note: op.basis field stores an instance but basis() returns the type B

# do-block syntax
function FunctionOperator(action::F, name::Symbol, basis::B; adjoint_action::Union{Function, Nothing}=nothing) where {F<:Function,B<:AbstractBasis}
    FunctionOperator{B}(name, basis, action, adjoint_action)
end

# Regular constructor
function FunctionOperator(name::Symbol, basis::B, action::F; adjoint_action::Union{Function, Nothing}=nothing) where {F<:Function,B<:AbstractBasis}
    FunctionOperator{B}(name, basis, action, adjoint_action)
end

# Apply
function Base.:*(op::FunctionOperator{B}, ket::Ket{B}) where {B<:AbstractBasis}
    op.action(ket)
end

# ProductKet lives in CompositeBasis{B1,B2}
function Base.:*(op::FunctionOperator{CompositeBasis{B1,B2}}, ket::ProductKet{B1,B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    op.action(ket)
end

# SumKet lives in CompositeBasis{B1,B2}
function Base.:*(op::FunctionOperator{CompositeBasis{B1,B2}}, ket::SumKet{B1,B2,T}) where {B1<:AbstractBasis, B2<:AbstractBasis, T}
    total = nothing
    for (k, w) in zip(ket.kets, ket.weights)
        result = op.action(k)
        # Check if result is zero BEFORE multiplying by weight (to avoid 0*symbolic issues)
        iszero(result) && continue
        term = result isa Number ? result * w : w * result
        total = isnothing(total) ? term : total + term
    end
    isnothing(total) ? 0 : total
end

function Base.:*(op::FunctionOperator, ket::WeightedKet)
    result = op * ket.ket
    result isa Number ? result * ket.weight : ket.weight * result
end

function Base.:*(op::FunctionOperator, ket::SumKet{B,T}) where {B,T}
    results = [op * Ket{B}(k.index) for k in ket.kets]
    weights = ket.weights
    
    total = nothing
    for (r, w) in zip(results, weights)
        term = r isa Number ? r * w : w * r
        if !iszero(term)
            total = isnothing(total) ? term : total + term
        end
    end
    isnothing(total) ? 0 : total
end

# Make callable
(op::FunctionOperator)(ket::AbstractKet) = op * ket

# Display
Base.show(io::IO, op::FunctionOperator) = print(io, op.name)

# Adjoint (wrapper)
@doc """
    AdjointFunctionOperator{B}

The Hermitian conjugate of a FunctionOperator. Created via `op'`.
""" AdjointFunctionOperator
struct AdjointFunctionOperator{B<:AbstractBasis} <: AbstractOperator{B}
    parent::FunctionOperator{B}
end

Base.adjoint(op::FunctionOperator) = AdjointFunctionOperator(op)
Base.adjoint(op::AdjointFunctionOperator) = op.parent

# basis inherited from AbstractOperator

Base.show(io::IO, op::AdjointFunctionOperator) = print(io, op.parent.name, "†")

# AdjointFunctionOperator application - uses adjoint_action if defined, otherwise symbolic
function Base.:*(op::AdjointFunctionOperator{B}, ket::Ket{B}) where B
    if !isnothing(op.parent.adjoint_action)
        op.parent.adjoint_action(ket)
    else
        OpKet(op, ket)
    end
end

# ProductKet lives in CompositeBasis{B1,B2}
function Base.:*(op::AdjointFunctionOperator{CompositeBasis{B1,B2}}, ket::ProductKet{B1,B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    if !isnothing(op.parent.adjoint_action)
        op.parent.adjoint_action(ket)
    else
        OpKet(op, ket)
    end
end

# SumKet lives in CompositeBasis{B1,B2}
function Base.:*(op::AdjointFunctionOperator{CompositeBasis{B1,B2}}, ket::SumKet{B1,B2,T}) where {B1<:AbstractBasis, B2<:AbstractBasis, T}
    total = nothing
    for (k, w) in zip(ket.kets, ket.weights)
        result = op * k
        term = result isa Number ? result * w : w * result
        if !iszero(term)
            total = isnothing(total) ? term : total + term
        end
    end
    isnothing(total) ? 0 : total
end

function Base.:*(op::AdjointFunctionOperator{B}, ket::WeightedKet{B}) where B
    result = op * ket.ket
    result isa Number ? result * ket.weight : ket.weight * result
end

function Base.:*(op::AdjointFunctionOperator{B}, ket::SumKet{B,T}) where {B,T}
    total = nothing
    for (k, w) in zip(ket.kets, ket.weights)
        result = op * k
        term = result isa Number ? result * w : w * result
        if !iszero(term)
            total = isnothing(total) ? term : total + term
        end
    end
    isnothing(total) ? 0 : total
end

# ============== Cross-type operations ==============

# FunctionOperator * Operator products
function Base.:*(op1::FunctionOperator{B}, op2::Operator{B}) where B
    OperatorProduct{B}([op1, op2])
end

function Base.:*(op1::Operator{B}, op2::FunctionOperator{B}) where B
    OperatorProduct{B}([op1, op2])
end

# AdjointFunctionOperator cross-type products
function Base.:*(op1::AdjointFunctionOperator{B}, op2::AbstractOperator{B}) where B
    OperatorProduct{B}([op1, op2])
end

function Base.:*(op1::AbstractOperator{B}, op2::AdjointFunctionOperator{B}) where B
    OperatorProduct{B}([op1, op2])
end

# Avoid ambiguity with OperatorProduct
function Base.:*(op1::AdjointFunctionOperator{B}, op2::OperatorProduct{B}) where B
    OperatorProduct{B}(vcat([op1], op2.ops))
end

function Base.:*(op1::OperatorProduct{B}, op2::AdjointFunctionOperator{B}) where B
    OperatorProduct{B}(vcat(op1.ops, [op2]))
end

# ============== Convenience constructors ==============

@doc """
    create_ladder_operators(basis; name=:a)

Create bosonic annihilation and creation operators for a Fock space basis.
Returns `(â, â†)` where:
- â|n⟩ = √n |n-1⟩  (annihilation)
- â†|n⟩ = √(n+1) |n+1⟩  (creation)

# Example
```julia
F = FockSpace(:mode)
Fb = Basis(F, :n)
â, â† = create_ladder_operators(Fb)

n0 = Ket(Fb, 0)
n1 = Ket(Fb, 1)

â * n1   # → √1 |0⟩
â† * n0  # → √1 |1⟩
```
""" create_ladder_operators
function create_ladder_operators(basis::B; name::Symbol=:a) where B<:AbstractBasis
    â = FunctionOperator(name, basis; 
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

# ============== Symbolic representations ==============

@doc """
    OpBra{B1, B2}

Symbolic representation of a bra with an operator applied: ⟨ψ|Ô
""" OpBra
struct OpBra{B1<:AbstractBasis, B2<:AbstractBasis}
    bra::AbstractBra{B1}
    op::AbstractOperator{B2}
end

Base.show(io::IO, ob::OpBra) = print(io, ob.bra, ob.op)

# Same-basis Bra * Operator
function Base.:*(bra::AbstractBra{B}, op::AbstractOperator{B}) where B
    OpBra(bra, op)
end

# Cross-basis Bra * Operator (same space, different bases)
function Base.:*(bra::AbstractBra{B1}, op::AbstractOperator{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    space(B1) == space(B2) || 
        throw(DimensionMismatch("Bra in $B1 and operator in $B2 are in different spaces"))
    OpBra(bra, op)
end

# ⟨ψ|Ô|ϕ⟩ = ⟨ψ|(Ô|ϕ⟩) - same basis
function Base.:*(ob::OpBra{B,B}, ket::AbstractKet{B}) where B
    result = ob.op * ket
    result isa Number ? result * (ob.bra * Ket{B}(nothing)) : ob.bra * result
end

# ⟨ψ|Ô|ϕ⟩ - cross-basis: apply operator to ket, then inner product
function Base.:*(ob::OpBra{B1,B2}, ket::AbstractKet{B3}) where {B1,B2,B3}
    # First apply operator to ket (may need basis transform)
    result = ob.op * ket
    # Then compute inner product with bra
    result isa Number ? result : ob.bra * result
end

# ⟨ψ|Ô¹Ô² = ⟨ψ|(Ô¹Ô²) - operator chaining
function Base.:*(ob::OpBra, op2::AbstractOperator)
    OpBra(ob.bra, ob.op * op2)
end

@doc """
    OpKet{B1, B2}

Symbolic representation of an operator applied to a ket: Ô|ψ⟩
Used when the result cannot be simplified.
""" OpKet
struct OpKet{B1<:AbstractBasis, B2<:AbstractBasis}
    op::AbstractOperator{B1}
    ket::AbstractKet{B2}
end

Base.show(io::IO, ok::OpKet) = print(io, ok.op, ok.ket)

# Cross-basis application helper - returns symbolic OpKet
function _cross_basis_apply(op::AbstractOperator{B1}, ket::AbstractKet{B2}) where {B1, B2}
    space(B1) == space(B2) || 
        throw(DimensionMismatch("Operator in $(B1) cannot act on ket in $(B2): different spaces"))
    OpKet(op, ket)
end

# ============== Identity operator ==============

@doc """
    IdentityOp{B}

The identity operator on basis B. Created via `IdentityOp(basis)`.
""" IdentityOp
struct IdentityOp{B<:AbstractBasis} <: AbstractOperator{B}
    function IdentityOp(::Type{B}) where B<:AbstractBasis
        new{B}()
    end
    function IdentityOp(b::B) where B<:AbstractBasis
        new{B}()
    end
end

# Constructor with explicit type parameter
IdentityOp{B}() where {B<:AbstractBasis} = IdentityOp(B)

# basis inherited from AbstractOperator

Base.:*(::IdentityOp{B}, ket::AbstractKet{B}) where B = ket
Base.:*(op::AbstractOperator{B}, ::IdentityOp{B}) where B = op
Base.:*(::IdentityOp{B}, op::AbstractOperator{B}) where B = op
Base.:*(::IdentityOp{B}, ::IdentityOp{B}) where B = IdentityOp{B}()
Base.:*(op::AdjointFunctionOperator{B}, ::IdentityOp{B}) where B = op

Base.adjoint(id::IdentityOp) = id

Base.show(io::IO, ::IdentityOp{B}) where B = print(io, "𝕀")

# ============== Zero checking ==============

Base.iszero(x::Number) = !(x isa AbstractSymbolic) && x == 0
Base.iszero(::AbstractKet) = false
Base.iszero(::AbstractOperator) = false

# ============== Tensor Product of Operators ==============

@doc """
    TensorOperator{B1,B2}(op1, op2)
    op1 ⊗ op2

Tensor product of two operators: Â ⊗ B̂. Acts on composite states as:
(Â ⊗ B̂)(|ψ⟩ ⊗ |ϕ⟩) = (Â|ψ⟩) ⊗ (B̂|ϕ⟩)

# Examples
```julia
H1, H2 = HilbertSpace(:A, 2), HilbertSpace(:B, 2)
B1, B2 = Basis(H1, :z), Basis(H2, :z)

# Create operators on each subsystem
σz_1 = up1 * up1' - down1 * down1'
σx_2 = up2 * down2' + down2 * up2'

# Tensor product
op = σz_1 ⊗ σx_2  # Acts on composite space

# Apply to product state
ket = up1 ⊗ up2
op * ket  # → (σz|↑⟩) ⊗ (σx|↑⟩)
```

See also: [`Operator`](@ref), [`IdentityOp`](@ref)
""" TensorOperator
struct TensorOperator{B1<:AbstractBasis, B2<:AbstractBasis, O1<:AbstractOperator{B1}, O2<:AbstractOperator{B2}} <: AbstractOperator{CompositeBasis{B1,B2}}
    op1::O1
    op2::O2
end

# Constructors
TensorOperator(op1::O1, op2::O2) where {B1, B2, O1<:AbstractOperator{B1}, O2<:AbstractOperator{B2}} = 
    TensorOperator{B1, B2, O1, O2}(op1, op2)

# Tensor product of operators
⊗(op1::AbstractOperator{B1}, op2::AbstractOperator{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis} = 
    TensorOperator(op1, op2)

# Display
function Base.show(io::IO, op::TensorOperator)
    print(io, "(", op.op1, ")⊗(", op.op2, ")")
end

# Apply TensorOperator to ProductKet: (Â⊗B̂)(|ψ⟩⊗|ϕ⟩) = (Â|ψ⟩)⊗(B̂|ϕ⟩)
function Base.:*(op::TensorOperator{B1,B2}, ket::ProductKet{B1,B2}) where {B1,B2}
    result1 = op.op1 * ket.ket1
    result2 = op.op2 * ket.ket2
    
    # Handle zero results
    (iszero(result1) || iszero(result2)) && return 0
    
    # Combine results based on their types
    _tensor_combine(result1, result2)
end

# Helper to combine tensor product results
function _tensor_combine(r1::Ket{B1}, r2::Ket{B2}) where {B1,B2}
    ProductKet(r1, r2)
end

function _tensor_combine(r1::WeightedKet{B1}, r2::Ket{B2}) where {B1,B2}
    r1.weight * ProductKet(r1.ket, r2)
end

function _tensor_combine(r1::Ket{B1}, r2::WeightedKet{B2}) where {B1,B2}
    r2.weight * ProductKet(r1, r2.ket)
end

function _tensor_combine(r1::WeightedKet{B1}, r2::WeightedKet{B2}) where {B1,B2}
    (r1.weight * r2.weight) * ProductKet(r1.ket, r2.ket)
end

function _tensor_combine(r1::SumKet{B1,T1}, r2::Ket{B2}) where {B1,B2,T1}
    kets = [ProductKet(k, r2) for k in r1.kets]
    SumKet(kets, r1.weights)
end

function _tensor_combine(r1::Ket{B1}, r2::SumKet{B2,T2}) where {B1,B2,T2}
    kets = [ProductKet(r1, k) for k in r2.kets]
    SumKet(kets, r2.weights)
end

function _tensor_combine(r1::WeightedKet{B1}, r2::SumKet{B2,T2}) where {B1,B2,T2}
    kets = [ProductKet(r1.ket, k) for k in r2.kets]
    SumKet(kets, r1.weight .* r2.weights)
end

function _tensor_combine(r1::SumKet{B1,T1}, r2::WeightedKet{B2}) where {B1,B2,T1}
    kets = [ProductKet(k, r2.ket) for k in r1.kets]
    SumKet(kets, r2.weight .* r1.weights)
end

function _tensor_combine(r1::SumKet{B1,T1}, r2::SumKet{B2,T2}) where {B1,B2,T1,T2}
    kets = ProductKet{B1,B2}[]
    weights = promote_type(T1,T2)[]
    for (k1, w1) in zip(r1.kets, r1.weights)
        for (k2, w2) in zip(r2.kets, r2.weights)
            push!(kets, ProductKet(k1, k2))
            push!(weights, w1 * w2)
        end
    end
    SumKet(kets, weights)
end

# Handle Number results from operator application
function _tensor_combine(r1::Number, r2)
    iszero(r1) ? 0 : r1 * _to_sum_ket(r2)
end

function _tensor_combine(r1, r2::Number)
    iszero(r2) ? 0 : r2 * _to_sum_ket(r1)
end

function _tensor_combine(r1::Number, r2::Number)
    r1 * r2
end

# Helper to convert various ket types to a consistent form
_to_sum_ket(k::Ket) = SumKet(k)
_to_sum_ket(k::WeightedKet) = SumKet(k)
_to_sum_ket(k::SumKet) = k

# Apply TensorOperator to SumKet
function Base.:*(op::TensorOperator{B1,B2}, ket::SumKet{B1,B2,T}) where {B1,B2,T}
    total = nothing
    for (k, w) in zip(ket.kets, ket.weights)
        result = op * k
        iszero(result) && continue
        term = result isa Number ? result * w : w * result
        total = isnothing(total) ? term : total + term
    end
    isnothing(total) ? 0 : total
end

# Adjoint of TensorOperator: (Â⊗B̂)† = Â†⊗B̂†
function Base.adjoint(op::TensorOperator{B1,B2}) where {B1,B2}
    TensorOperator(adjoint(op.op1), adjoint(op.op2))
end

# Scalar multiplication
function Base.:*(a::Number, op::TensorOperator{B1,B2}) where {B1,B2}
    ScaledOperator{CompositeBasis{B1,B2}, typeof(a)}(a, op)
end

# Addition of TensorOperators
function Base.:+(op1::TensorOperator{B1,B2}, op2::TensorOperator{B1,B2}) where {B1,B2}
    SumOperator{CompositeBasis{B1,B2}}([op1, op2])
end

function Base.:+(op1::TensorOperator{B1,B2}, op2::AbstractOperator{CompositeBasis{B1,B2}}) where {B1,B2}
    terms2 = op2 isa SumOperator ? op2.terms : [op2]
    SumOperator{CompositeBasis{B1,B2}}(vcat([op1], terms2))
end

function Base.:+(op1::AbstractOperator{CompositeBasis{B1,B2}}, op2::TensorOperator{B1,B2}) where {B1,B2}
    terms1 = op1 isa SumOperator ? op1.terms : [op1]
    SumOperator{CompositeBasis{B1,B2}}(vcat(terms1, [op2]))
end

# Multiplication of TensorOperators: (Â⊗B̂)(Ĉ⊗D̂) = (ÂĈ)⊗(B̂D̂)
function Base.:*(op1::TensorOperator{B1,B2}, op2::TensorOperator{B1,B2}) where {B1,B2}
    TensorOperator(op1.op1 * op2.op1, op1.op2 * op2.op2)
end

# TensorOperator with other composite operators
function Base.:*(op1::TensorOperator{B1,B2}, op2::AbstractOperator{CompositeBasis{B1,B2}}) where {B1,B2}
    OperatorProduct{CompositeBasis{B1,B2}}([op1, op2])
end

function Base.:*(op1::AbstractOperator{CompositeBasis{B1,B2}}, op2::TensorOperator{B1,B2}) where {B1,B2}
    OperatorProduct{CompositeBasis{B1,B2}}([op1, op2])
end

# Identity tensor products: 𝕀⊗Â and Â⊗𝕀
⊗(::IdentityOp{B1}, op2::AbstractOperator{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis} = TensorOperator(IdentityOp{B1}(), op2)
⊗(op1::AbstractOperator{B1}, ::IdentityOp{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis} = TensorOperator(op1, IdentityOp{B2}())
⊗(::IdentityOp{B1}, ::IdentityOp{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis} = TensorOperator(IdentityOp{B1}(), IdentityOp{B2}())

# IdentityOp on ProductKet
function Base.:*(::IdentityOp{CompositeBasis{B1,B2}}, ket::ProductKet{B1,B2}) where {B1,B2}
    ket
end

function Base.:*(::IdentityOp{CompositeBasis{B1,B2}}, ket::SumKet{B1,B2,T}) where {B1,B2,T}
    ket
end

# TensorOperator with IdentityOp simplification
function Base.:*(op::TensorOperator{B1,B2,IdentityOp{B1},O2}, ket::ProductKet{B1,B2}) where {B1,B2,O2}
    result2 = op.op2 * ket.ket2
    iszero(result2) && return 0
    _tensor_combine(ket.ket1, result2)
end

function Base.:*(op::TensorOperator{B1,B2,O1,IdentityOp{B2}}, ket::ProductKet{B1,B2}) where {B1,B2,O1}
    result1 = op.op1 * ket.ket1
    iszero(result1) && return 0
    _tensor_combine(result1, ket.ket2)
end

# ============== Lifted Operators (single-system operator acting on composite) ==============

@doc """
    lift(op::AbstractOperator{B}, position::Int, target_basis::CompositeBasis)
    lift(op, :first, B1 ⊗ B2)   # Â → Â ⊗ 𝕀
    lift(op, :second, B1 ⊗ B2)  # B̂ → 𝕀 ⊗ B̂

Lift a single-system operator to act on a composite system by tensoring with identity.

# Examples
```julia
H1, H2 = HilbertSpace(:A, 2), HilbertSpace(:B, 2)
B1, B2 = Basis(H1, :z), Basis(H2, :z)

σz = up * up' - down * down'  # Operator on B1

# Lift to composite: σz ⊗ 𝕀
σz_lifted = lift(σz, :first, B1 ⊗ B2)
# Or equivalently:
σz_lifted = lift(σz, 1, CompositeBasis{typeof(B1), typeof(B2)})
```
""" lift
function lift(op::AbstractOperator{B1}, position::Int, ::Type{CompositeBasis{BA,BB}}) where {B1,BA,BB}
    if position == 1
        B1 == BA || throw(ArgumentError("Operator basis $B1 doesn't match first component $BA"))
        return TensorOperator(op, IdentityOp{BB}())
    elseif position == 2
        B1 == BB || throw(ArgumentError("Operator basis $B1 doesn't match second component $BB"))
        return TensorOperator(IdentityOp{BA}(), op)
    else
        throw(ArgumentError("Position must be 1 or 2 for two-component composite basis"))
    end
end

lift(op::AbstractOperator{B}, position::Symbol, cb::CompositeBasis{BA,BB}) where {B,BA,BB} = 
    lift(op, position == :first ? 1 : 2, CompositeBasis{BA,BB})

lift(op::AbstractOperator{B}, position::Int, cb::CompositeBasis{BA,BB}) where {B,BA,BB} = 
    lift(op, position, CompositeBasis{BA,BB})

# ============== Swap / Reorder Operations ==============

@doc """
    swap(ket::ProductKet)
    swap(bra::ProductBra)
    swap(op::TensorOperator)

Swap the order of subsystems in a tensor product.
|ψ⟩⊗|ϕ⟩ → |ϕ⟩⊗|ψ⟩
Â⊗B̂ → B̂⊗Â

# Examples
```julia
ket = up1 ⊗ down2
swap(ket)  # → down2 ⊗ up1

op = σz ⊗ σx
swap(op)   # → σx ⊗ σz
```
""" swap
function swap(ket::ProductKet{B1,B2}) where {B1,B2}
    ProductKet(ket.ket2, ket.ket1)
end

function swap(bra::ProductBra{B1,B2}) where {B1,B2}
    ProductBra(bra.bra2, bra.bra1)
end

function swap(sk::SumKet{B1,B2,T}) where {B1,B2,T}
    swapped_kets = [ProductKet(k.ket2, k.ket1) for k in sk.kets]
    SumKet(swapped_kets, sk.weights; name=sk.display_name)
end

function swap(sb::SumBra{B1,B2,T}) where {B1,B2,T}
    swapped_bras = [ProductBra(b.bra2, b.bra1) for b in sb.bras]
    SumBra(swapped_bras, sb.weights; name=sb.display_name)
end

function swap(op::TensorOperator{B1,B2}) where {B1,B2}
    TensorOperator(op.op2, op.op1)
end

@doc """
    reorder(obj, target_basis::CompositeBasis)
    reorder(obj, basis1 ⊗ basis2)

Reorder a tensor product (ket, bra, or operator) to match the given target basis order.
If the object's basis order matches, returns as-is. If swapped, applies swap.

# Examples
```julia
H1, H2 = HilbertSpace(:A, 2), HilbertSpace(:B, 2)
B1, B2 = Basis(H1, :z), Basis(H2, :z)

ket = Ket(B2, :↑) ⊗ Ket(B1, :↓)  # In B2⊗B1 order
target = B1 ⊗ B2

reorder(ket, target)  # → |↓⟩⊗|↑⟩ in B1⊗B2 order

op = σx_B2 ⊗ σz_B1  # Operator in B2⊗B1
reorder(op, target)  # → σz_B1 ⊗ σx_B2 in B1⊗B2 order
```
""" reorder
function reorder(ket::ProductKet{A1,A2}, ::Type{CompositeBasis{B1,B2}}) where {A1,A2,B1,B2}
    if A1 == B1 && A2 == B2
        return ket
    elseif A1 == B2 && A2 == B1
        return swap(ket)
    else
        throw(ArgumentError("Cannot reorder: bases don't match. Got $A1⊗$A2, target $B1⊗$B2"))
    end
end

function reorder(sk::SumKet{A1,A2,T}, ::Type{CompositeBasis{B1,B2}}) where {A1,A2,T,B1,B2}
    if A1 == B1 && A2 == B2
        return sk
    elseif A1 == B2 && A2 == B1
        return swap(sk)
    else
        throw(ArgumentError("Cannot reorder: bases don't match. Got $A1⊗$A2, target $B1⊗$B2"))
    end
end

function reorder(bra::ProductBra{A1,A2}, ::Type{CompositeBasis{B1,B2}}) where {A1,A2,B1,B2}
    if A1 == B1 && A2 == B2
        return bra
    elseif A1 == B2 && A2 == B1
        return swap(bra)
    else
        throw(ArgumentError("Cannot reorder: bases don't match. Got $A1⊗$A2, target $B1⊗$B2"))
    end
end

function reorder(sb::SumBra{A1,A2,T}, ::Type{CompositeBasis{B1,B2}}) where {A1,A2,T,B1,B2}
    if A1 == B1 && A2 == B2
        return sb
    elseif A1 == B2 && A2 == B1
        return swap(sb)
    else
        throw(ArgumentError("Cannot reorder: bases don't match. Got $A1⊗$A2, target $B1⊗$B2"))
    end
end

function reorder(op::TensorOperator{A1,A2}, ::Type{CompositeBasis{B1,B2}}) where {A1,A2,B1,B2}
    if A1 == B1 && A2 == B2
        return op
    elseif A1 == B2 && A2 == B1
        return swap(op)
    else
        throw(ArgumentError("Cannot reorder: bases don't match. Got $A1⊗$A2, target $B1⊗$B2"))
    end
end

# Convenience: reorder to a CompositeBasis instance
reorder(obj, cb::CompositeBasis{B1,B2}) where {B1,B2} = reorder(obj, CompositeBasis{B1,B2})

# ============== Partial Trace ==============

@doc """
    partial_trace(op::TensorOperator, subsystem::Int)

Compute the partial trace over one subsystem of a tensor product operator.
Returns an operator on the remaining subsystem.

- `partial_trace(Â⊗B̂, 1)` traces over first subsystem → Tr(Â)·B̂
- `partial_trace(Â⊗B̂, 2)` traces over second subsystem → Tr(B̂)·Â

Note: For general operators, the trace is symbolic unless the operator is 
in a form where trace can be computed (e.g., projectors |ψ⟩⟨ψ|).

# Examples
```julia
# For projectors: Tr(|ψ⟩⟨ψ|) = 1
P1 = up * up'
P2 = down * down'
op = P1 ⊗ P2

partial_trace(op, 1)  # → Tr(P1)·P2 = P2
partial_trace(op, 2)  # → Tr(P2)·P1 = P1
```
""" partial_trace
function partial_trace(op::TensorOperator{B1,B2,Operator{B1,T1},O2}, subsystem::Int) where {B1,B2,T1,O2}
    if subsystem == 1
        # Trace over first: Tr(|ψ⟩⟨ϕ|) = ⟨ϕ|ψ⟩
        inner = op.op1.bra * Ket{B1}(op.op1.ket.index)
        coeff = op.op1.coeff * inner
        return coeff * op.op2
    elseif subsystem == 2
        return partial_trace_second(op)
    else
        throw(ArgumentError("subsystem must be 1 or 2"))
    end
end

function partial_trace(op::TensorOperator{B1,B2,O1,Operator{B2,T2}}, subsystem::Int) where {B1,B2,O1,T2}
    if subsystem == 2
        # Trace over second: Tr(|ψ⟩⟨ϕ|) = ⟨ϕ|ψ⟩
        inner = op.op2.bra * Ket{B2}(op.op2.ket.index)
        coeff = op.op2.coeff * inner
        return coeff * op.op1
    elseif subsystem == 1
        return partial_trace_first(op)
    else
        throw(ArgumentError("subsystem must be 1 or 2"))
    end
end

# Both are Operator type
function partial_trace(op::TensorOperator{B1,B2,Operator{B1,T1},Operator{B2,T2}}, subsystem::Int) where {B1,B2,T1,T2}
    if subsystem == 1
        inner = op.op1.bra * Ket{B1}(op.op1.ket.index)
        coeff = op.op1.coeff * inner
        return coeff * op.op2
    elseif subsystem == 2
        inner = op.op2.bra * Ket{B2}(op.op2.ket.index)
        coeff = op.op2.coeff * inner
        return coeff * op.op1
    else
        throw(ArgumentError("subsystem must be 1 or 2"))
    end
end

# Fallback for general operators - returns symbolic
function partial_trace(op::TensorOperator{B1,B2}, subsystem::Int) where {B1,B2}
    throw(ArgumentError("Partial trace not implemented for general operator types. Use Operator (outer product) form."))
end

# ============== Automatic basis transformation ==============
#
# When an operator in basis B1 acts on a ket in basis B2,
# and B1 ≠ B2 but they're in the same space, automatically
# transform the ket to B1 before applying the operator.
#

"""
    _try_transform_and_apply(op, ket) -> result or nothing

Try to transform ket to operator's basis and apply. Returns nothing if no transform exists.
"""
function _try_transform_and_apply(op::AbstractOperator{B1}, ket::AbstractKet{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    # Check if same space
    space(B1) == space(B2) || return nothing
    # Check if ket is already in operator's basis (for composite bases)
    # ProductKet{A1,A2} is in CompositeBasis{A1,A2}, SumKet{A1,A2} is also in CompositeBasis{A1,A2}
    ket_basis = _ket_composite_basis(ket)
    if !isnothing(ket_basis) && ket_basis == B1
        return nothing  # Already in same basis - let specific same-basis method handle it
    end
    # Check if transform exists
    has_transform(B2, B1) || return nothing
    # Transform ket to operator's basis and apply
    transformed_ket = transform(ket, B1)
    return op * transformed_ket
end

# Helper to get the composite basis type from a product ket
_ket_composite_basis(::ProductKet{A1,A2}) where {A1,A2} = CompositeBasis{A1,A2}
_ket_composite_basis(::SumKet{A1,A2}) where {A1,A2} = CompositeBasis{A1,A2}
_ket_composite_basis(::AbstractKet) = nothing

# Cross-basis application for all operator types
# These catch cases where operator basis B1 ≠ ket basis B2

# General AbstractOperator × AbstractKet (fallback with auto-transform)
function Base.:*(op::AbstractOperator{B1}, ket::AbstractKet{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    result = _try_transform_and_apply(op, ket)
    if isnothing(result)
        if space(B1) == space(B2)
            throw(ArgumentError("No transform registered from $B2 to $B1. Define one with define_transform!"))
        else
            throw(ArgumentError("Cannot apply operator on $(space(B1)) to ket on $(space(B2))"))
        end
    end
    return result
end
