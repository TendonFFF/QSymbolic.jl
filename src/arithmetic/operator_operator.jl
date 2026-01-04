# Operator * Operator multiplication

# No additional exports - uses base * operator

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

# ============== Identity Operator Methods ==============

# Apply identity
Base.:*(::Identity, ket::AbstractKet) = ket

# Operator algebra with identity
Base.:*(::Identity{S}, op::AbstractOperator{S}) where S = op
Base.:*(op::AbstractOperator{S}, ::Identity{S}) where S = op
Base.:*(id1::Identity{S}, ::Identity{S}) where S = id1

# Adjoint
Base.adjoint(id::Identity) = id

# Scalar multiplication of Identity creates OperatorSum (c * 𝕀)
Base.:*(c::Number, id::Identity{S}) where S = OperatorSum([id], [c])
Base.:*(id::Identity{S}, c::Number) where S = c * id

# ============== OperatorSum Arithmetic ==============

# Addition: OperatorSum + anything
Base.:+(opsum::OperatorSum{S}, op::AbstractOperator{S}) where S = 
    OperatorSum(vcat(opsum.operators, [op]), vcat(opsum.weights, [1]))

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
