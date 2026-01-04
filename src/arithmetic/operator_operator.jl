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

# ============== FunctionOperator * Outer/Operator ==============
# F * |a⟩⟨b| = (F|a⟩)⟨b| - apply F to the ket part

function Base.:*(fop::FunctionOperator{S}, outer::Outer{S}) where S
    # F * |a⟩⟨b| = (F|a⟩) ⊗ ⟨b|
    new_ket = fop * outer.ket
    if iszero(new_ket)
        return 0
    end
    # new_ket could be Ket, WeightedKet, or SumKet
    return new_ket * outer.bra
end

function Base.:*(outer::Outer{S}, fop::FunctionOperator{S}) where S
    # |a⟩⟨b| * F = |a⟩ ⊗ (F†|b⟩)†
    # Need the adjoint of F to act on the bra side
    # ⟨b|F = (F†|b⟩)†
    fop_adj = fop'  # This may throw if adjoint not defined
    ket_from_bra = Ket(outer.bra.index)  # Convert bra index back to ket
    # Actually we need to construct a proper Ket in the same basis as the bra
    B = typeof(outer.bra).parameters[1]  # Extract basis type
    ket_equiv = Ket(B(), outer.bra.index)
    new_ket = fop_adj * ket_equiv
    if iszero(new_ket)
        return 0
    end
    # Result is |a⟩⟨new_bra| where new_bra = new_ket'
    return outer.ket * new_ket'
end

# FunctionOperator * Operator (Operator is a sum of Outers)
function Base.:*(fop::FunctionOperator{S}, op::Operator{S}) where S
    results = []
    for (outer, w) in zip(op.outers, op.weights)
        result = fop * outer
        if !iszero(result)
            push!(results, w * result)
        end
    end
    isempty(results) && return 0
    # Sum all results
    total = results[1]
    for i in 2:length(results)
        total = total + results[i]
    end
    total
end

function Base.:*(op::Operator{S}, fop::FunctionOperator{S}) where S
    results = []
    for (outer, w) in zip(op.outers, op.weights)
        result = outer * fop
        if !iszero(result)
            push!(results, w * result)
        end
    end
    isempty(results) && return 0
    # Sum all results
    total = results[1]
    for i in 2:length(results)
        total = total + results[i]
    end
    total
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

# Generic fallback: any two AbstractOperators in same space → OperatorSum
# More specific methods (OperatorSum+X, X+Identity) take precedence
Base.:+(op1::AbstractOperator{S}, op2::AbstractOperator{S}) where {S<:AbstractSpace} = OperatorSum([op1, op2], [1, 1])

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

# ============== OperatorSum * Operator/Outer (Same Space - Evaluate Eagerly) ==============

# Helper to check if result is zero (handles operators which can't be zero)
_is_zero_result(x) = iszero(x)
_is_zero_result(::AbstractOperator) = false  # Operators are never "zero" in this sense

# Helper to check if a result is an OperatorProduct (indicating no contraction occurred)
_is_lazy_product(::OperatorProduct) = true
_is_lazy_product(_) = false

# OperatorSum * Outer: distribute and contract if possible
function Base.:*(opsum::OperatorSum{S}, outer::Outer{S}) where S
    results = []
    weights_out = []
    
    for (op, w) in zip(opsum.operators, opsum.weights)
        result = op * outer
        # If result is OperatorProduct, no contraction happened - fall back to lazy
        if _is_lazy_product(result)
            return OperatorProduct([opsum, outer])
        end
        if !_is_zero_result(result)
            push!(results, result)
            push!(weights_out, w)
        end
    end
    
    isempty(results) && return 0
    
    # Sum all the contracted results
    total = weights_out[1] * results[1]
    for i in 2:length(results)
        total = total + weights_out[i] * results[i]
    end
    total
end

function Base.:*(outer::Outer{S}, opsum::OperatorSum{S}) where S
    results = []
    weights_out = []
    
    for (op, w) in zip(opsum.operators, opsum.weights)
        result = outer * op
        # If result is OperatorProduct, no contraction happened - fall back to lazy
        if _is_lazy_product(result)
            return OperatorProduct([outer, opsum])
        end
        if !_is_zero_result(result)
            push!(results, result)
            push!(weights_out, w)
        end
    end
    
    isempty(results) && return 0
    
    # Sum all the contracted results
    total = weights_out[1] * results[1]
    for i in 2:length(results)
        total = total + weights_out[i] * results[i]
    end
    total
end

# OperatorSum * Operator: similar treatment
function Base.:*(opsum::OperatorSum{S}, op2::Operator{S}) where S
    results = []
    weights_out = []
    
    for (op, w) in zip(opsum.operators, opsum.weights)
        result = op * op2
        if _is_lazy_product(result)
            return OperatorProduct([opsum, op2])
        end
        if !_is_zero_result(result)
            push!(results, result)
            push!(weights_out, w)
        end
    end
    
    isempty(results) && return 0
    
    total = weights_out[1] * results[1]
    for i in 2:length(results)
        total = total + weights_out[i] * results[i]
    end
    total
end

function Base.:*(op1::Operator{S}, opsum::OperatorSum{S}) where S
    results = []
    weights_out = []
    
    for (op, w) in zip(opsum.operators, opsum.weights)
        result = op1 * op
        if _is_lazy_product(result)
            return OperatorProduct([op1, opsum])
        end
        if !_is_zero_result(result)
            push!(results, result)
            push!(weights_out, w)
        end
    end
    
    isempty(results) && return 0
    
    total = weights_out[1] * results[1]
    for i in 2:length(results)
        total = total + weights_out[i] * results[i]
    end
    total
end

# ============== Generic Lazy Multiplication (OperatorProduct) ==============
# For all other operator combinations, create a lazy OperatorProduct

# Generic fallback for same space: any two AbstractOperators → OperatorProduct
# More specific methods (Outer*Outer, Operator*Operator, etc.) take precedence
Base.:*(op1::AbstractOperator{S}, op2::AbstractOperator{S}) where S = OperatorProduct([op1, op2])

# Generic fallback for different spaces: create lazy OperatorProduct
# Validity is checked when applied to a ket
Base.:*(op1::AbstractOperator, op2::AbstractOperator) = OperatorProduct([op1, op2])

# Flatten nested OperatorProducts (no space constraint)
Base.:*(op1::OperatorProduct, op2::AbstractOperator) = 
    OperatorProduct(vcat(op1.operators, [op2]))
Base.:*(op1::AbstractOperator, op2::OperatorProduct) = 
    OperatorProduct(vcat([op1], op2.operators))
Base.:*(op1::OperatorProduct, op2::OperatorProduct) = 
    OperatorProduct(vcat(op1.operators, op2.operators))

# ============== OperatorProduct Application ==============
# Apply operators right-to-left: (ABC)|ψ⟩ = A(B(C|ψ⟩))

function Base.:*(opprod::OperatorProduct, ket::AbstractKet)
    result = ket
    for op in reverse(opprod.operators)
        result = op * result
    end
    result
end
