# Operator types: Outer, Operator, Identity, FunctionOperator

# Quantum operators - Restructured to work with bra-ket arithmetic

# space() - get the space an operator acts on
space(::AbstractOperator{S}) where S = S

# ============== Outer Product |ψ⟩⟨ϕ| ==============

@doc """
    Outer{S<:AbstractSpace}(ket, bra)
    ket * bra'  (automatically creates Outer)

Basic outer product operator |ψ⟩⟨ϕ|. Allows cross-basis operators
where the ket and bra can be in different bases, as long as they're
in the same space.

Application uses bra-ket contraction:
    (|ψ⟩⟨ϕ|)|χ⟩ = ⟨ϕ|χ⟩ |ψ⟩

# Examples
```julia
H, Hb = HilbertSpace(:H, 2)
up = Ket(Hb, :↑)
down = Ket(Hb, :↓)

# Projector onto |↑⟩
P_up = up * up'  # |↑⟩⟨↑|

# Apply: P_up|↑⟩ = |↑⟩, P_up|↓⟩ = 0
P_up * up   # → |↑⟩
P_up * down # → 0

# Ladder operator |↑⟩⟨↓|
σ_plus = up * down'
σ_plus * down  # → |↑⟩
```

See also: [`Operator`](@ref), [`Identity`](@ref), [`FunctionOperator`](@ref)
"""
struct Outer{S<:AbstractSpace} <: AbstractOperator{S}
    ket::AbstractKet
    bra::AbstractBra
    space::S
    
    function Outer(ket::AbstractKet, bra::AbstractBra)
        # Verify same space
        ket_space = space(basis(ket))
        bra_space = space(basis(bra))
        ket_space == bra_space || throw(ArgumentError("Ket and bra must be in the same space"))
        new{typeof(ket_space)}(ket, bra, ket_space)
    end
end

# Create outer product from ket * bra'
Base.:*(ket::AbstractKet, bra::AbstractBra) = Outer(ket, bra)

# Display
function Base.show(io::IO, op::Outer)
    print(io, _ket_str(op.ket), _bra_str(op.bra))
end

# Helper functions for display
_ket_str(k::Ket) = "|$(k.index)⟩"
_ket_str(k::WeightedKet) = "$(k.weight)·|$(k.ket.index)⟩"
_ket_str(k::ProductKet) = "|" * join([string(ki.index) for ki in k.kets], "⊗") * "⟩"
_ket_str(k::SumKet) = "(" * join(["$(w)·|$(ki.index)⟩" for (ki, w) in zip(k.kets, k.weights)], " + ") * ")"
_ket_str(k::AbstractKet) = "|ket⟩"

_bra_str(b::Bra) = "⟨$(b.index)|"
_bra_str(b::WeightedBra) = "$(b.weight)·⟨$(b.bra.index)|"
_bra_str(b::ProductBra) = "⟨" * join([string(bi.index) for bi in b.bras], "⊗") * "|"
_bra_str(b::SumBra) = "(" * join(["$(w)·⟨$(bi.index)|" for (bi, w) in zip(b.bras, b.weights)], " + ") * ")"