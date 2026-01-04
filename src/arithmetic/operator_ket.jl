# Operator * Ket applications

# No additional exports - uses base * operator

# Apply outer product to basic Ket: (|ψ⟩⟨ϕ|)|χ⟩ = ⟨ϕ|χ⟩ |ψ⟩
function Base.:*(op::Outer, ket::Ket)
    # Use bra-ket contraction
    inner = op.bra * ket
    if !(inner isa AbstractSymbolic) && iszero(inner)
        return 0
    else
        return inner * op.ket
    end
end

# Apply outer product to WeightedKet
function Base.:*(op::Outer, wk::WeightedKet)
    result = op * wk.ket
    if result isa Number && iszero(result)
        return 0
    end
    return wk.weight * result
end

# Apply outer product to SumKet
function Base.:*(op::Outer, sk::SumKet)
    total = nothing
    for (ket, w) in zip(sk.kets, sk.weights)
        result = op * ket
        if result isa Number && iszero(result)
            continue
        end
        term = w * result
        total = isnothing(total) ? term : total + term
    end
    isnothing(total) ? 0 : total
end

# Adjoint: (|ψ⟩⟨ϕ|)† = |ϕ⟩⟨ψ|
Base.adjoint(op::Outer) = Outer(adjoint(op.bra), adjoint(op.ket))

# Apply operator to basic Ket: decomposes to bra-ket contractions
function Base.:*(op::Operator, ket::Ket)
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

# Apply Operator to WeightedKet
function Base.:*(op::Operator, wk::WeightedKet)
    result = op * wk.ket
    if result isa Number && iszero(result)
        return 0
    end
    return wk.weight * result
end

# Apply Operator to SumKet
function Base.:*(op::Operator, sk::SumKet)
    total = nothing
    for (ket, w) in zip(sk.kets, sk.weights)
        result = op * ket
        if result isa Number && iszero(result)
            continue
        end
        term = w * result
        total = isnothing(total) ? term : total + term
    end
    isnothing(total) ? 0 : total
end

# Apply FunctionOperator to basic Ket
function Base.:*(op::FunctionOperator{S, B}, ket::Ket{B}) where {S, B}
    op.action(ket)
end

# Apply FunctionOperator to WeightedKet: preserve weight
function Base.:*(op::FunctionOperator{S, B}, wk::WeightedKet{B}) where {S, B}
    result = op.action(wk.ket)
    if result isa Number && iszero(result)
        return 0
    end
    return wk.weight * result
end

# Apply FunctionOperator to SumKet: distribute over terms
function Base.:*(op::FunctionOperator{S, B}, sk::SumKet{B}) where {S, B}
    total = nothing
    for (ket, w) in zip(sk.kets, sk.weights)
        result = op.action(ket)
        if result isa Number && iszero(result)
            continue
        end
        term = w * result
        total = isnothing(total) ? term : total + term
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

# Scalar multiplication for Operator
function Base.:*(a::Number, op::Operator{S}) where S
    Operator{S}(op.outers, [a * w for w in op.weights])
end

function Base.:*(a::Number, outer::Outer{S}) where S
    Operator{S}([outer], [a])
end

Base.:*(op::AbstractOperator, a::Number) = a * op
Base.:/(op::AbstractOperator, a::Number) = (1/a) * op

# ==================== PARTIAL APPLICATION ON PRODUCT STATES ====================
# Operators acting on one subsystem of a composite (product) state

# Helper: Find the index of the ket in ProductKet that matches the operator's space
function _find_matching_subsystem(pk::ProductKet, ::Type{S}) where S<:AbstractSpace
    for (i, k) in enumerate(pk.kets)
        ket_space = space(basis(typeof(k)))
        if ket_space == S
            return i
        end
    end
    return nothing
end

# Apply Outer to ProductKet: (|ψ⟩⟨ϕ|)_{A} (|a⟩⊗|b⟩) = ⟨ϕ|a⟩ |ψ⟩⊗|b⟩  (if A is subsystem of |a⟩)
function Base.:*(op::Outer{S}, pk::ProductKet) where S
    idx = _find_matching_subsystem(pk, S)
    if isnothing(idx)
        throw(ArgumentError("Operator space $S not found in ProductKet subsystems"))
    end
    
    # Extract the matching ket and the rest
    target_ket = pk.kets[idx]
    other_kets = [pk.kets[i] for i in 1:length(pk.kets) if i != idx]
    
    # Apply outer product: ⟨ϕ|target⟩ |ψ⟩
    inner = op.bra * target_ket
    
    if !(inner isa AbstractSymbolic) && iszero(inner)
        return 0
    end
    
    # Result: inner * |ψ⟩ ⊗ (other kets)
    result_ket = op.ket
    
    if isempty(other_kets)
        return inner * result_ket
    else
        # Reconstruct product state
        new_product = length(other_kets) == 1 ? 
            result_ket ⊗ other_kets[1] :
            result_ket ⊗ ProductKet(other_kets)
        return inner * new_product
    end
end

# Apply Operator to ProductKet (sum of Outers)
function Base.:*(op::Operator{S}, pk::ProductKet) where S
    idx = _find_matching_subsystem(pk, S)
    if isnothing(idx)
        throw(ArgumentError("Operator space $S not found in ProductKet subsystems"))
    end
    
    total = nothing
    for (outer, w) in zip(op.outers, op.weights)
        result = outer * pk
        if !(result isa AbstractSymbolic) && iszero(result)
            continue
        end
        term = w * result
        total = isnothing(total) ? term : total + term
    end
    isnothing(total) ? 0 : total
end

# Apply FunctionOperator to ProductKet
function Base.:*(op::FunctionOperator{S, B}, pk::ProductKet) where {S, B}
    idx = _find_matching_subsystem(pk, S)
    if isnothing(idx)
        throw(ArgumentError("Operator space $S not found in ProductKet subsystems"))
    end
    
    # Extract the matching ket and the rest
    target_ket = pk.kets[idx]
    other_kets = [pk.kets[i] for i in 1:length(pk.kets) if i != idx]
    
    # Apply function operator to matching ket
    result = op.action(target_ket)
    
    if result isa Number && iszero(result)
        return 0
    end
    
    if isempty(other_kets)
        return result
    else
        # Result could be a Ket, WeightedKet, or SumKet
        # Tensor with remaining kets
        if result isa Ket
            return length(other_kets) == 1 ?
                result ⊗ other_kets[1] :
                result ⊗ ProductKet(other_kets)
        elseif result isa WeightedKet
            new_product = length(other_kets) == 1 ?
                result.ket ⊗ other_kets[1] :
                result.ket ⊗ ProductKet(other_kets)
            return result.weight * new_product
        else
            # For SumKet or other types, return as-is (may need extension)
            return result ⊗ (length(other_kets) == 1 ? other_kets[1] : ProductKet(other_kets))
        end
    end
end

# Apply FunctionOperator to WeightedKet with CompositeBasis (e.g., from Outer * ProductKet)
function Base.:*(op::FunctionOperator{S, B}, wk::WeightedKet{CB}) where {S, B, CB<:CompositeBasis}
    # The underlying ket should be a ProductKet - apply op to it, then rescale
    result = op * wk.ket
    if result isa Number && iszero(result)
        return 0
    end
    return wk.weight * result
end

# Apply FunctionOperator to SumKet with CompositeBasis
function Base.:*(op::FunctionOperator{S, B}, sk::SumKet{CB}) where {S, B, CB<:CompositeBasis}
    total = nothing
    for (ket, w) in zip(sk.kets, sk.weights)
        result = op * ket
        if result isa Number && iszero(result)
            continue
        end
        term = w * result
        total = isnothing(total) ? term : total + term
    end
    isnothing(total) ? 0 : total
end
