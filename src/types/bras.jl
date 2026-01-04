# Bra types: Bra, WeightedBra, SumBra, ProductBra

# ==================== BRA TYPES ====================
# Bras are lazy adjoints of kets with the same structure

@doc """
    Bra{B<:AbstractBasis}(index)
    Bra(basis, index)

Basic bra state ⟨index| in a specific basis.
Usually created via adjoint: `ket'`.

**Note**: A basis is **required**. You must explicitly provide a basis object.

See also: [`Ket`](@ref), [`Basis`](@ref)
""" Bra
struct Bra{B<:AbstractBasis} <: AbstractBra{B}
    index::KetIndex

    function Bra(::B, idx::KetIndex=nothing) where B<:AbstractBasis
        new{B}(idx)
    end
    function Bra(::B, idx::Int) where B<:AbstractBasis
        new{B}(Symbol(idx))
    end
    function Bra(::B, idx::Tuple{Vararg{Int}}) where B<:AbstractBasis
        new{B}(Tuple(Symbol(i) for i in idx))
    end
    function Bra{B}(idx::Union{Symbol, Int, Nothing}=nothing) where B<:AbstractBasis
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function Bra{B}(idx::AbstractSymbolic) where B<:AbstractBasis
        new{B}(idx)
    end
    function Bra{B}(idx::Tuple{Vararg{SingleIndexValue}}) where B<:AbstractBasis
        new{B}(idx)
    end
end

@doc """
    ProductBra{Bs<:Tuple}(bras::Vector{Bra})
    ProductBra(bra1, bra2, ...)

Tensor product of bras: ⟨ψ₁|⊗⟨ψ₂|⊗...⊗⟨ψₙ|.
Created via adjoint of ProductKet.

**Note**: ProductBras are **order-independent** (bosonic/symmetric),
matching the behavior of ProductKet.

See also: [`ProductKet`](@ref), [`Bra`](@ref)
""" ProductBra
struct ProductBra{Bs<:Tuple} <: AbstractBra{CompositeBasis{Bs}}
    bras::Vector{Bra}
    
    function ProductBra(bras::Vector{Bra})
        if length(bras) < 2
            throw(ArgumentError("ProductBra requires at least 2 bras"))
        end
        # Canonicalize order: sort bras by basis type for order-independence
        sorted_bras = _canonicalize_product_bras(bras)
        basis_types = Tuple{[basis(typeof(b)) for b in sorted_bras]...}
        new{basis_types}(sorted_bras)
    end
    
    function ProductBra(bra1::Bra, bra2::Bra, bras::Bra...)
        ProductBra([bra1, bra2, bras...])
    end
end

# Canonicalize bra order for bosonic (order-independent) behavior
function _canonicalize_product_bras(bras::Vector{Bra})
    # Sort by basis type, then by index for deterministic ordering
    return sort(bras, by = b -> (string(basis(typeof(b))), string(b.index)))
end

@doc """
    WeightedBra{B,T}(bra, weight)

A bra multiplied by a scalar weight.
Created via adjoint of WeightedKet (with complex conjugate weight).

See also: [`WeightedKet`](@ref)
""" WeightedBra
struct WeightedBra{B<:AbstractBasis, T} <: AbstractBra{B}
    bra::Union{Bra{B}, ProductBra}
    weight::T
    
    function WeightedBra(bra::Bra{B}, weight::T) where {B<:AbstractBasis, T}
        new{B, T}(bra, weight)
    end
    function WeightedBra(pb::ProductBra{Bs}, weight::T) where {Bs,T}
        CB = CompositeBasis{Bs}
        new{CB, T}(pb, weight)
    end
end

@doc """
    SumBra{B,T}(bras, weights; name=nothing)

Linear combination of bras with weights.
Created via adjoint of SumKet (with complex conjugate weights).

See also: [`SumKet`](@ref)
""" SumBra
struct SumBra{B<:AbstractBasis, T} <: AbstractBra{B}
    display_name::Union{Symbol,Nothing}
    bras::Vector{Union{Bra{B}, ProductBra}}
    weights::Vector{T}

    function SumBra(bras::AbstractVector, weights::AbstractVector{T}; name::Union{Symbol, Nothing}=nothing) where T
        if isempty(bras)
            throw(ArgumentError("Cannot create SumBra with empty bra vector"))
        end
        B = basis(typeof(bras[1]))
        for b in bras
            basis(typeof(b)) == B || throw(ArgumentError("All bras must be in the same basis"))
        end
        new{B, T}(name, bras, weights)
    end
    
    function SumBra(bra::Union{Bra{B}, ProductBra{Bs}}; name::Union{Symbol, Nothing}=nothing) where {B, Bs}
        if bra isa Bra
            new{B, Int}(name, [bra], [1])
        else  # ProductBra
            CB = CompositeBasis{Bs}
            new{CB, Int}(name, [bra], [1])
        end
    end
    
    function SumBra(wb::WeightedBra{B, T}; name::Union{Symbol, Nothing}=nothing) where {B, T}
        new{B, T}(name, [wb.bra], [wb.weight])
    end
end
