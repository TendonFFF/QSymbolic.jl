# Quantum operators with user-defined actions

"""
    Operator{B}(name, basis, action)

A quantum operator acting on states in basis `B`. Created via `Operator(name, basis) do ket ... end`.

The `action` function receives a `BasisKet` and should return:
- A number (scalar)
- A ket (`BasisKet`, `weightedKet`, or `sumKet`)
- Zero (for annihilation of vacuum, etc.)

# Example
```julia
H = HilbertSpace(:H, 2)
Zb = Basis(H, :z)

# Pauli Z operator: σz|↑⟩ = |↑⟩, σz|↓⟩ = -|↓⟩
σz = Operator(:σz, Zb) do ket
    ket.index == :↑ ? BasisKet(Zb, :↑) : -1 * BasisKet(Zb, :↓)
end

# Apply to a state
σz * BasisKet(Zb, :↑)  # → |↑⟩
```

See also: [`AdjointOperator`](@ref), [`OperatorProduct`](@ref)
"""
struct Operator{B<:AbstractBasis}
    name::Symbol
    basis::B
    action::Function
end

"""
    basis(op::Operator)

Get the basis that an operator is defined in.
"""
basis(op::Operator) = op.basis

# do-block syntax: Operator(name, basis) do ket ... end
Operator(action::Function, name::Symbol, basis::B) where B<:AbstractBasis = Operator{B}(name, basis, action)

# Apply operator to kets
"""
    op * ket

Apply operator `op` to ket. The operator's action function determines the result.
Operator and ket must be in the same basis.
"""
function Base.:*(op::Operator{B}, ket::BasisKet{B}) where {B<:AbstractBasis}
    op.action(ket)
end

# Cross-basis application: requires basis transform
function Base.:*(op::Operator{B1}, ket::BasisKet{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    # Check if spaces are compatible
    space(B1) == space(B2) || 
        throw(DimensionMismatch("Operator in $(B1) cannot act on ket in $(B2): different spaces"))
    # For now, return symbolic OpKet for cross-basis application
    # TODO: auto-transform if transform is defined
    OpKet(op, ket)
end

function Base.:*(op::Operator, ket::weightedKet)
    result = op * ket.Ket
    result isa Number ? result * ket.weight : ket.weight * result
end

function Base.:*(op::Operator, ket::sumKet{B,T}) where {B,T}
    results = [op * BasisKet{B}(k.index) for k in ket.kets]
    weights = ket.weights
    
    # Combine results, handling numbers and kets
    total = nothing
    for (r, w) in zip(results, weights)
        term = r isa Number ? r * w : w * r
        total = isnothing(total) ? term : total + term
    end
    total
end

# Make operator callable
(op::Operator)(ket::AbstractKet) = op * ket

# Display
Base.show(io::IO, op::Operator) = print(io, op.name)

# Adjoint operator (wrapper)
"""
    AdjointOperator{B}

The Hermitian conjugate of an operator. Created via `op'`.
"""
struct AdjointOperator{B<:AbstractBasis}
    parent::Operator{B}
end

Base.adjoint(op::Operator) = AdjointOperator(op)
Base.adjoint(op::AdjointOperator) = op.parent

basis(op::AdjointOperator) = basis(op.parent)

Base.show(io::IO, op::AdjointOperator) = print(io, op.parent.name, "†")

# Apply adjoint operator to bras (from the right)
function Base.:*(bra::AbstractBra, op::Operator)
    # ⟨ψ|Ô returns a new bra-like object
    # For now, return symbolic
    OpBra(bra, op)
end

"""
    OpBra

Symbolic representation of a bra with an operator applied: ⟨ψ|Ô
"""
struct OpBra{B1<:AbstractBasis, B2<:AbstractBasis}
    bra::AbstractBra{B1}
    op::Operator{B2}
end

Base.show(io::IO, ob::OpBra) = print(io, ob.bra, ob.op.name)

# ⟨ψ|Ô|ϕ⟩ = ⟨ψ|(Ô|ϕ⟩)
function Base.:*(ob::OpBra, ket::AbstractKet)
    result = ob.op * ket
    result isa Number ? result * (ob.bra * BasisKet{basis(ob.bra)}(nothing)) : ob.bra * result
end

"""
    OpKet

Symbolic representation of an operator applied to a ket: Ô|ψ⟩
Used when the result cannot be simplified (e.g., cross-basis application).
"""
struct OpKet{B1<:AbstractBasis, B2<:AbstractBasis}
    op::Operator{B1}
    ket::AbstractKet{B2}
end

Base.show(io::IO, ok::OpKet) = print(io, ok.op.name, ok.ket)

# Operator algebra (symbolic for now)
"""
    OperatorProduct{B}

Product of operators in the same basis: ÂB̂
"""
struct OperatorProduct{B<:AbstractBasis}
    ops::Vector{Union{Operator{B}, AdjointOperator{B}}}
end

basis(prod::OperatorProduct{B}) where B = B

function Base.:*(op1::Operator{B}, op2::Operator{B}) where B
    OperatorProduct{B}([op1, op2])
end

function Base.:*(prod::OperatorProduct{B}, op::Operator{B}) where B
    OperatorProduct{B}([prod.ops..., op])
end

function Base.:*(op::Operator{B}, prod::OperatorProduct{B}) where B
    OperatorProduct{B}([op, prod.ops...])
end

function Base.show(io::IO, prod::OperatorProduct)
    for op in prod.ops
        print(io, op)
    end
end

# Apply product of operators: (ÂB̂)|ψ⟩ = Â(B̂|ψ⟩)
function Base.:*(prod::OperatorProduct{B}, ket::BasisKet{B}) where B
    result = ket
    for op in reverse(prod.ops)
        result = op * result
    end
    result
end

function Base.:*(prod::OperatorProduct, ket::AbstractKet)
    result = ket
    for op in reverse(prod.ops)
        result = op * result
    end
    result
end
