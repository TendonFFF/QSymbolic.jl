# OperatorSum container for lazy operator addition

# OperatorSum container for lazy operator addition

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
