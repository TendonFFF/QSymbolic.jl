# Operator * Ket applications

# No additional exports - uses base * operator

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
"""
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
