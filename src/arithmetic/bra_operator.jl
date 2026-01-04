# Bra * Operator contractions

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
