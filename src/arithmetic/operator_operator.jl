# Operator * Operator multiplication

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

# Identity operator multiplication
Base.:*(::Identity{S}, op::AbstractOperator{S}) where S = op
Base.:*(op::AbstractOperator{S}, ::Identity{S}) where S = op

# ============== Identity Operator ==============

@doc """
    Identity{S<:AbstractSpace}(space)

The identity operator on space S. Basis-independent.

# Examples
```julia
H, Hb = HilbertSpace(:H, 2)
I = Identity(H)

ψ = Ket(Hb, :ψ)
I * ψ  # → |ψ⟩
```

See also: [`Operator`](@ref), [`Outer`](@ref)
"""
struct Identity{S<:AbstractSpace} <: AbstractOperator{S}
    space::S
    
    function Identity(s::S) where S<:AbstractSpace
        new{S}(s)
    end
end

# Apply identity
Base.:*(::Identity, ket::AbstractKet) = ket

# Operator algebra with identity
Base.:*(::Identity{S}, op::AbstractOperator{S}) where S = op
Base.:*(op::AbstractOperator{S}, ::Identity{S}) where S = op
Base.:*(::Identity{S}, ::Identity{S}) where S = Identity(S)

# Adjoint
Base.adjoint(id::Identity) = id

# Display
Base.show(io::IO, ::Identity) = print(io, "𝕀")

# ============== Function-based Operator ==============

@doc """
    FunctionOperator{S<:AbstractSpace, B<:AbstractBasis}(basis, action; adjoint_action=nothing, name=:F)
    FunctionOperator(basis) do ket ... end

Function-based operator that applies a user-defined function to kets in a specific basis.
The function receives a ket in the operator's basis and returns any AbstractKet.

If the input ket is in a different basis, the operator automatically applies a basis
transform before applying the action (if a transform is defined).

Constructed with do-block syntax:
```julia
F_op = FunctionOperator(basis) do ket
    # Process ket in 'basis', return AbstractKet or Number
    ...
end
```

Optionally provide `adjoint_action` for the adjoint operator and `name` for display.

# Examples
```julia
F, Fb = HilbertSpace(:Fock, Inf)

# Annihilation operator: â|n⟩ = √n |n-1⟩
â = FunctionOperator(Fb) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end

# With adjoint (creation operator):
â = FunctionOperator(Fb; 
    adjoint_action = ket -> begin
        n = parse(Int, string(ket.index))
        √(n + 1) * Ket(Fb, n + 1)
    end,
    name = :â
) do ket
    n = parse(Int, string(ket.index))
    n == 0 ? 0 : √n * Ket(Fb, n - 1)
end

# Cross-basis: if ket is in different basis, transform is applied first
Fb2 = Basis(F, :energy)
# Define transform between bases first
define_transform!(typeof(Fb2), typeof(Fb)) do k
    # ... transformation logic ...
end
â * Ket(Fb2, :E0)  # Transforms to Fb basis first, then applies â
```

See also: [`Operator`](@ref), [`Outer`](@ref), [`define_transform!`](@ref)
"""
struct FunctionOperator{S<:AbstractSpace, B<:AbstractBasis} <: AbstractOperator{S}
    space::S
    basis::B
    action::Function
    adjoint_action::Union{Function, Nothing}
    name::Symbol
    
    function FunctionOperator{S,B}(space::S, basis::B, action::Function, adjoint_action::Union{Function, Nothing}, name::Symbol) where {S<:AbstractSpace, B<:AbstractBasis}
        space(basis) == space || throw(ArgumentError("Basis must be in the same space as operator"))
        new{S,B}(space, basis, action, adjoint_action, name)
    end
end

# Constructor with do-block
function FunctionOperator(action::F, basis::B; adjoint_action::Union{Function, Nothing}=nothing, name::Symbol=:F) where {F<:Function, B<:AbstractBasis}
    s = space(basis)
    FunctionOperator{typeof(s),typeof(basis)}(s, basis, action, adjoint_action, name)
end

# Regular constructor
function FunctionOperator(basis::B, action::F; adjoint_action::Union{Function, Nothing}=nothing, name::Symbol=:F) where {F<:Function, B<:AbstractBasis}
    s = space(basis)
    FunctionOperator{typeof(s),typeof(basis)}(s, basis, action, adjoint_action, name)
end
