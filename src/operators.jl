# Quantum operators

# ============== Abstract Operator Type ==============

"""
    AbstractOperator{B}

Abstract supertype for all quantum operators in basis `B`.
"""
abstract type AbstractOperator{B<:AbstractBasis} end

# ============== Outer Product Operator |ψ⟩⟨ϕ| ==============

"""
    Operator{B}(ket, bra)
    ket * bra'  (automatically creates Operator)

An operator formed as the outer product |ψ⟩⟨ϕ|. This is the standard
quantum mechanical representation of operators.

When applied to a state: (|ψ⟩⟨ϕ|)|χ⟩ = |ψ⟩⟨ϕ|χ⟩ = ⟨ϕ|χ⟩ |ψ⟩

# Examples
```julia
H = HilbertSpace(:H, 2)
Zb = Basis(H, :z)
up = BasisKet(Zb, :↑)
down = BasisKet(Zb, :↓)

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
"""
struct Operator{B<:AbstractBasis, T} <: AbstractOperator{B}
    ket::AbstractKet{B}
    bra::AbstractBra{B}
    coeff::T
    
    function Operator(ket::AbstractKet{B}, bra::AbstractBra{B}, coeff::T=1) where {B<:AbstractBasis, T}
        new{B, T}(ket, bra, coeff)
    end
end

"""
    basis(op::AbstractOperator)

Get the basis type that an operator is defined in.
"""
basis(::AbstractOperator{B}) where B = B

"""
    space(op::AbstractOperator)

Get the space type that an operator acts on.
"""
space(::AbstractOperator{B}) where B = space(B)

# Create operator from ket * bra
function Base.:*(ket::BasisKet{B}, bra::BasisBra{B}) where B<:AbstractBasis
    Operator(ket, bra, 1)
end

function Base.:*(ket::weightedKet{B}, bra::BasisBra{B}) where B
    Operator(ket.Ket, bra, ket.weight)
end

function Base.:*(ket::BasisKet{B}, bra::weightedBra{B}) where B
    Operator(ket, bra.Bra, bra.weight)
end

function Base.:*(ket::weightedKet{B}, bra::weightedBra{B}) where B
    Operator(ket.Ket, bra.Bra, ket.weight * bra.weight)
end

# Display
function Base.show(io::IO, op::Operator)
    if op.coeff == 1
        print(io, _ket_str(op.ket), _bra_str(op.bra))
    else
        print(io, op.coeff, "·", _ket_str(op.ket), _bra_str(op.bra))
    end
end

# Helper functions for display
_ket_str(k::BasisKet) = "|$(k.index)⟩"
_bra_str(b::BasisBra) = "⟨$(b.index)|"

# Apply operator to ket: (|ψ⟩⟨ϕ|)|χ⟩ = ⟨ϕ|χ⟩ |ψ⟩
function Base.:*(op::Operator{B}, ket::BasisKet{B}) where B
    inner = BasisBra{B}(op.bra.index) * ket
    if inner == 0
        return 0
    else
        op.coeff * inner * op.ket
    end
end

function Base.:*(op::Operator{B}, ket::weightedKet{B}) where B
    result = op * ket.Ket
    result isa Number ? result * ket.weight : ket.weight * result
end

function Base.:*(op::Operator{B}, ket::sumKet{B,T}) where {B,T}
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
    new_ket = BasisKet{B}(op.bra.index)
    new_bra = BasisBra{B}(op.ket isa BasisKet ? op.ket.index : nothing)
    Operator(new_ket, new_bra, conj(op.coeff))
end

# ============== Sum of Operators ==============

"""
    SumOperator{B}

Sum of operators: Â + B̂. Created via operator addition.
"""
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
function Base.:*(op::SumOperator{B}, ket::BasisKet{B}) where B
    total = nothing
    for term in op.terms
        result = term * ket
        if !iszero(result)
            total = isnothing(total) ? result : total + result
        end
    end
    isnothing(total) ? 0 : total
end

function Base.:*(op::SumOperator{B}, ket::weightedKet{B}) where B
    total = nothing
    for term in op.terms
        result = term * ket
        if !iszero(result)
            total = isnothing(total) ? result : total + result
        end
    end
    isnothing(total) ? 0 : total
end

function Base.:*(op::SumOperator{B}, ket::sumKet{B,T}) where {B,T}
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

# SumOperator on SumProductKet (CompositeBasis)
function Base.:*(op::SumOperator{CompositeBasis{B1,B2}}, ket::SumProductKet{B1,B2,T}) where {B1,B2,T}
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

"""
    ScaledOperator{B,T}

Operator multiplied by a scalar. Created via `scalar * operator`.
"""
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

function Base.:*(op::ScaledOperator{B}, ket::BasisKet{B}) where B
    result = op.op * ket
    result isa Number ? op.coeff * result : op.coeff * result
end

function Base.:*(op::ScaledOperator{B}, ket::weightedKet{B}) where B
    result = op.op * ket
    result isa Number ? op.coeff * result : op.coeff * result
end

function Base.:*(op::ScaledOperator{B}, ket::sumKet{B,T}) where {B,T}
    result = op.op * ket
    result isa Number ? op.coeff * result : op.coeff * result
end

# ScaledOperator on ProductKet (CompositeBasis)
function Base.:*(op::ScaledOperator{CompositeBasis{B1,B2}}, ket::ProductKet{B1,B2}) where {B1,B2}
    result = op.op * ket
    result isa Number ? op.coeff * result : op.coeff * result
end

# ScaledOperator on SumProductKet (CompositeBasis)
function Base.:*(op::ScaledOperator{CompositeBasis{B1,B2}}, ket::SumProductKet{B1,B2,T}) where {B1,B2,T}
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

"""
    OperatorProduct{B}

Product of operators: ÂB̂. Created via operator multiplication.
"""
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
function Base.:*(prod::OperatorProduct{B}, ket::BasisKet{B}) where B
    result = ket
    for op in reverse(prod.ops)
        result = op * result
        iszero(result) && return 0
    end
    result
end

function Base.:*(prod::OperatorProduct{B}, ket::weightedKet{B}) where B
    result = ket
    for op in reverse(prod.ops)
        result = op * result
        iszero(result) && return 0
    end
    result
end

function Base.:*(prod::OperatorProduct{B}, ket::sumKet{B,T}) where {B,T}
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

# OperatorProduct on SumProductKet (CompositeBasis)
function Base.:*(prod::OperatorProduct{CompositeBasis{B1,B2}}, ket::SumProductKet{B1,B2,T}) where {B1,B2,T}
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

"""
    FunctionOperator{B}(name, basis, action; adjoint_action=nothing)

A quantum operator defined by a function. Created via `FunctionOperator(name, basis) do ket ... end`.

The `action` function receives a `BasisKet` and should return:
- A number (scalar)
- A ket (`BasisKet`, `weightedKet`, or `sumKet`)
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
    n == 0 ? 0 : √n * BasisKet(Fb, n - 1)
end

# With explicit adjoint action (creation operator):
â = FunctionOperator(:â, Fb; 
    adjoint_action = ket -> begin
        n = parse(Int, string(ket.index))
        √(n + 1) * BasisKet(Fb, n + 1)
    end
) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * BasisKet(Fb, n - 1)
end
```

See also: [`Operator`](@ref), [`create_ladder_operators`](@ref)
"""
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
function Base.:*(op::FunctionOperator{B}, ket::BasisKet{B}) where {B<:AbstractBasis}
    op.action(ket)
end

# ProductKet lives in CompositeBasis{B1,B2}
function Base.:*(op::FunctionOperator{CompositeBasis{B1,B2}}, ket::ProductKet{B1,B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    op.action(ket)
end

# SumProductKet lives in CompositeBasis{B1,B2}
function Base.:*(op::FunctionOperator{CompositeBasis{B1,B2}}, ket::SumProductKet{B1,B2,T}) where {B1<:AbstractBasis, B2<:AbstractBasis, T}
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

function Base.:*(op::FunctionOperator, ket::weightedKet)
    result = op * ket.Ket
    result isa Number ? result * ket.weight : ket.weight * result
end

function Base.:*(op::FunctionOperator, ket::sumKet{B,T}) where {B,T}
    results = [op * BasisKet{B}(k.index) for k in ket.kets]
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
"""
    AdjointFunctionOperator{B}

The Hermitian conjugate of a FunctionOperator. Created via `op'`.
"""
struct AdjointFunctionOperator{B<:AbstractBasis} <: AbstractOperator{B}
    parent::FunctionOperator{B}
end

Base.adjoint(op::FunctionOperator) = AdjointFunctionOperator(op)
Base.adjoint(op::AdjointFunctionOperator) = op.parent

# basis inherited from AbstractOperator

Base.show(io::IO, op::AdjointFunctionOperator) = print(io, op.parent.name, "†")

# AdjointFunctionOperator application - uses adjoint_action if defined, otherwise symbolic
function Base.:*(op::AdjointFunctionOperator{B}, ket::BasisKet{B}) where B
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

# SumProductKet lives in CompositeBasis{B1,B2}
function Base.:*(op::AdjointFunctionOperator{CompositeBasis{B1,B2}}, ket::SumProductKet{B1,B2,T}) where {B1<:AbstractBasis, B2<:AbstractBasis, T}
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

function Base.:*(op::AdjointFunctionOperator{B}, ket::weightedKet{B}) where B
    result = op * ket.Ket
    result isa Number ? result * ket.weight : ket.weight * result
end

function Base.:*(op::AdjointFunctionOperator{B}, ket::sumKet{B,T}) where {B,T}
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

"""
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

n0 = BasisKet(Fb, 0)
n1 = BasisKet(Fb, 1)

â * n1   # → √1 |0⟩
â† * n0  # → √1 |1⟩
```
"""
function create_ladder_operators(basis::B; name::Symbol=:a) where B<:AbstractBasis
    â = FunctionOperator(name, basis; 
        adjoint_action = ket -> begin
            n = ket.index isa Symbol ? parse(Int, string(ket.index)) : Int(ket.index)
            √(n + 1) * BasisKet(basis, n + 1)
        end
    ) do ket
        n = ket.index isa Symbol ? parse(Int, string(ket.index)) : Int(ket.index)
        n == 0 ? 0 : √n * BasisKet(basis, n - 1)
    end
    
    return â, â'
end

# ============== Symbolic representations ==============

"""
    OpBra{B1, B2}

Symbolic representation of a bra with an operator applied: ⟨ψ|Ô
"""
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
    result isa Number ? result * (ob.bra * BasisKet{B}(nothing)) : ob.bra * result
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

"""
    OpKet{B1, B2}

Symbolic representation of an operator applied to a ket: Ô|ψ⟩
Used when the result cannot be simplified.
"""
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

"""
    IdentityOp{B}

The identity operator on basis B. Created via `IdentityOp(basis)`.
"""
struct IdentityOp{B<:AbstractBasis} <: AbstractOperator{B}
    function IdentityOp(::Type{B}) where B<:AbstractBasis
        new{B}()
    end
    function IdentityOp(b::B) where B<:AbstractBasis
        new{B}()
    end
end

# basis inherited from AbstractOperator

Base.:*(::IdentityOp{B}, ket::AbstractKet{B}) where B = ket
Base.:*(op::AbstractOperator{B}, ::IdentityOp{B}) where B = op
Base.:*(::IdentityOp{B}, op::AbstractOperator{B}) where B = op
Base.:*(::IdentityOp{B}, ::IdentityOp{B}) where B = IdentityOp{B}()

Base.adjoint(id::IdentityOp) = id

Base.show(io::IO, ::IdentityOp{B}) where B = print(io, "𝕀")

# ============== Zero checking ==============

Base.iszero(x::Number) = x == 0
Base.iszero(::AbstractKet) = false
Base.iszero(::AbstractOperator) = false

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
    # ProductKet{A1,A2} is in CompositeBasis{A1,A2}, SumProductKet{A1,A2} is also in CompositeBasis{A1,A2}
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
_ket_composite_basis(::SumProductKet{A1,A2}) where {A1,A2} = CompositeBasis{A1,A2}
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
