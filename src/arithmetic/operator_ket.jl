# Operator * Ket applications

# No additional exports - uses base * operator

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

# Apply operator to ket: decomposes to bra-ket contractions
function Base.:*(op::Operator, ket::AbstractKet)
    total = nothing
    for (outer, w) in zip(op.outers, op.weights)
        result = outer * ket
        if !(result isa AbstractSymbolic) && iszero(result)
            continue
        end
        term = w * result
        total = isnothing(total) ? term : total + term
    end
    isnothing(total) ? 0 : total
end

# Apply FunctionOperator to ket
function Base.:*(op::FunctionOperator{S, B}, ket::AbstractKet{B}) where {S, B}
    op.action(ket)
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

# Scalar multiplication for Operator
function Base.:*(a::Number, op::Operator{S}) where S
    Operator{S}(op.outers, [a * w for w in op.weights])
end

function Base.:*(a::Number, outer::Outer{S}) where S
    Operator{S}([outer], [a])
end

Base.:*(op::AbstractOperator, a::Number) = a * op
Base.:/(op::AbstractOperator, a::Number) = (1/a) * op
