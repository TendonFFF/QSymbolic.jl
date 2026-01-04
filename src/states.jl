# New Ket and Bra state vectors with improved hierarchy
# This is a complete restructuring of the state type system

# ==================== KET TYPES ====================

# Index types - support arbitrary number of indices
# An index can be:
# - A single value: Symbol, Nothing, or AbstractSymbolic
# - A tuple of any length: (n,), (n, m), (n, m, k), etc.
const SingleIndexValue = Union{Symbol, Nothing, AbstractSymbolic}
const KetIndex = Union{SingleIndexValue, Tuple{Vararg{SingleIndexValue}}}

@doc """
    Ket{B<:AbstractBasis}(index)
    Ket(basis, index)
    Ket(space, index)  # uses DefaultBasis

Basic ket state |index⟩ in a specific basis. This is the fundamental building block.

Index can be:
- A `Symbol` (e.g., `:↑`, `:ψ`)
- An `Int` (converted to Symbol)
- A symbolic expression (`Sym`, `SymExpr`)
- A tuple for multi-index (e.g., `(:n, :m)` for composite bases)

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> Zb = Basis(H, :z);

julia> up = Ket(Zb, :↑)   # |↑⟩ in z-basis
|↑⟩

julia> ψ = Ket(H, :ψ)     # |ψ⟩ in default basis
|ψ⟩
```

See also: [`Bra`](@ref), [`WeightedKet`](@ref), [`SumKet`](@ref), [`ProductKet`](@ref)
""" Ket
struct Ket{B<:AbstractBasis} <: AbstractKet{B}
    index::KetIndex

    function Ket(::B, idx::KetIndex=nothing) where B<:AbstractBasis
        new{B}(idx)
    end
    function Ket(::B, idx::Int) where B<:AbstractBasis
        new{B}(Symbol(idx))
    end
    function Ket(::B, idx::Tuple{Vararg{Int}}) where B<:AbstractBasis
        # Convert all Int elements to Symbol
        new{B}(Tuple(Symbol(i) for i in idx))
    end
    function Ket{B}(idx::Union{Symbol, Int, Nothing}=nothing) where B<:AbstractBasis
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function Ket{B}(idx::AbstractSymbolic) where B<:AbstractBasis
        new{B}(idx)
    end
    function Ket{B}(idx::Tuple{Vararg{SingleIndexValue}}) where B<:AbstractBasis
        new{B}(idx)
    end
    function Ket(space::AbstractSpace, idx::Union{Symbol, Int, Nothing}=nothing)
        B = DefaultBasis{typeof(space)}
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function Ket(space::AbstractSpace, idx::AbstractSymbolic)
        B = DefaultBasis{typeof(space)}
        new{B}(idx)
    end
    function Ket(space::AbstractSpace, idx::Tuple{Vararg{SingleIndexValue}})
        B = DefaultBasis{typeof(space)}
        new{B}(idx)
    end
    function Ket(space::AbstractSpace, idx::Tuple{Vararg{Int}})
        B = DefaultBasis{typeof(space)}
        new{B}(Tuple(Symbol(i) for i in idx))
    end
end

@doc """
    ProductKet{B1<:AbstractBasis, B2<:AbstractBasis}(ket1, ket2)

Tensor product of two basic kets: |ψ⟩⊗|ϕ⟩.
Lives in `CompositeBasis{B1,B2}`. Same level as Ket in the hierarchy.

# Examples
```jldoctest
julia> H1, H2 = HilbertSpace(:A, 2), HilbertSpace(:B, 2);

julia> ψ, ϕ = Ket(H1, :ψ), Ket(H2, :ϕ);

julia> ψ ⊗ ϕ
|ψ⟩⊗|ϕ⟩
```

See also: [`Ket`](@ref), [`⊗`](@ref)
""" ProductKet
struct ProductKet{B1<:AbstractBasis, B2<:AbstractBasis} <: AbstractKet{CompositeBasis{B1,B2}}
    ket1::Ket{B1}
    ket2::Ket{B2}
end

@doc """
    WeightedKet{B<:AbstractBasis, T}(ket, weight)

A ket (Ket or ProductKet) multiplied by a scalar weight.
Created automatically via scalar multiplication: `weight * ket`.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ = Ket(H, :ψ);

julia> 2 * ψ
2·|ψ⟩
```

See also: [`Ket`](@ref), [`SumKet`](@ref)
""" WeightedKet
struct WeightedKet{B<:AbstractBasis, T} <: AbstractKet{B}
    ket::Union{Ket{B}, ProductKet}  # Can wrap either basic type
    weight::T
    
    # Constructor that handles both Ket and ProductKet
    function WeightedKet(ket::Ket{B}, weight::T) where {B<:AbstractBasis, T}
        new{B, T}(ket, weight)
    end
    function WeightedKet(pk::ProductKet{B1,B2}, weight::T) where {B1,B2,T}
        CB = CompositeBasis{B1,B2}
        new{CB, T}(pk, weight)
    end
end

@doc """
    SumKet{B<:AbstractBasis, T}(kets, weights; name=nothing)

Linear combination of kets with weights: Σᵢ wᵢ|ψᵢ⟩.
Created automatically via addition/subtraction of kets.

The kets can be either `Ket` or `ProductKet` (but all must share the same basis structure).

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ, ϕ = Ket(H, :ψ), Ket(H, :ϕ);

julia> ψ + ϕ
|ψ⟩ + |ϕ⟩

julia> (ψ + ϕ) / √2
0.7071067811865475·|ψ⟩ + 0.7071067811865475·|ϕ⟩
```

See also: [`Ket`](@ref), [`WeightedKet`](@ref)
""" SumKet
struct SumKet{B<:AbstractBasis, T} <: AbstractKet{B}
    display_name::Union{Symbol,Nothing}
    kets::Vector{Union{Ket{B}, ProductKet}}  # Vector of basic kets
    weights::Vector{T}

    function SumKet(kets::AbstractVector, weights::AbstractVector{T}; name::Union{Symbol, Nothing}=nothing) where T
        # Infer basis from first ket
        if isempty(kets)
            throw(ArgumentError("Cannot create SumKet with empty ket vector"))
        end
        B = basis(typeof(kets[1]))
        # Ensure all kets are in same basis
        for k in kets
            basis(typeof(k)) == B || throw(ArgumentError("All kets must be in the same basis"))
        end
        new{B, T}(name, kets, weights)
    end
    
    # Convenience constructor from single ket
    function SumKet(ket::Union{Ket{B}, ProductKet{B1,B2}}; name::Union{Symbol, Nothing}=nothing) where {B, B1, B2}
        if ket isa Ket
            new{B, Int}(name, [ket], [1])
        else  # ProductKet
            CB = CompositeBasis{B1, B2}
            new{CB, Int}(name, [ket], [1])
        end
    end
    
    # Constructor from WeightedKet
    function SumKet(wk::WeightedKet{B, T}; name::Union{Symbol, Nothing}=nothing) where {B, T}
        new{B, T}(name, [wk.ket], [wk.weight])
    end
end

# ==================== BRA TYPES ====================
# Bras are lazy adjoints of kets with the same structure

@doc """
    Bra{B<:AbstractBasis}(index)
    Bra(basis, index)
    Bra(space, index)  # uses DefaultBasis

Basic bra state ⟨index| in a specific basis.
Usually created via adjoint: `ket'`.

See also: [`Ket`](@ref)
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
    function Bra(space::AbstractSpace, idx::Union{Symbol, Int, Nothing}=nothing)
        B = DefaultBasis{typeof(space)}
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function Bra(space::AbstractSpace, idx::AbstractSymbolic)
        B = DefaultBasis{typeof(space)}
        new{B}(idx)
    end
    function Bra(space::AbstractSpace, idx::Tuple{Vararg{SingleIndexValue}})
        B = DefaultBasis{typeof(space)}
        new{B}(idx)
    end
    function Bra(space::AbstractSpace, idx::Tuple{Vararg{Int}})
        B = DefaultBasis{typeof(space)}
        new{B}(Tuple(Symbol(i) for i in idx))
    end
end

@doc """
    ProductBra{B1,B2}(bra1, bra2)

Tensor product of two bras: ⟨ψ|⊗⟨ϕ|.
Created via adjoint of ProductKet.

See also: [`ProductKet`](@ref)
""" ProductBra
struct ProductBra{B1<:AbstractBasis, B2<:AbstractBasis} <: AbstractBra{CompositeBasis{B1,B2}}
    bra1::Bra{B1}
    bra2::Bra{B2}
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
    function WeightedBra(pb::ProductBra{B1,B2}, weight::T) where {B1,B2,T}
        CB = CompositeBasis{B1,B2}
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
    
    function SumBra(bra::Union{Bra{B}, ProductBra{B1,B2}}; name::Union{Symbol, Nothing}=nothing) where {B, B1, B2}
        if bra isa Bra
            new{B, Int}(name, [bra], [1])
        else  # ProductBra
            CB = CompositeBasis{B1, B2}
            new{CB, Int}(name, [bra], [1])
        end
    end
    
    function SumBra(wb::WeightedBra{B, T}; name::Union{Symbol, Nothing}=nothing) where {B, T}
        new{B, T}(name, [wb.bra], [wb.weight])
    end
end

# ==================== UTILITY FUNCTIONS ====================

# basis() and space() functions
basis(::Ket{B}) where B = B
basis(::Bra{B}) where B = B
basis(::Type{Ket{B}}) where B = B
basis(::Type{Bra{B}}) where B = B
basis(::ProductKet{B1,B2}) where {B1,B2} = CompositeBasis{B1,B2}
basis(::ProductBra{B1,B2}) where {B1,B2} = CompositeBasis{B1,B2}
basis(::Type{ProductKet{B1,B2}}) where {B1,B2} = CompositeBasis{B1,B2}
basis(::Type{ProductBra{B1,B2}}) where {B1,B2} = CompositeBasis{B1,B2}
basis(::WeightedKet{B}) where B = B
basis(::WeightedBra{B}) where B = B
basis(::Type{WeightedKet{B,T}}) where {B,T} = B
basis(::Type{WeightedBra{B,T}}) where {B,T} = B
basis(::SumKet{B}) where B = B
basis(::SumBra{B}) where B = B
basis(::Type{SumKet{B,T}}) where {B,T} = B
basis(::Type{SumBra{B,T}}) where {B,T} = B

space(k::Union{Ket, Bra, ProductKet, ProductBra, WeightedKet, WeightedBra, SumKet, SumBra}) = space(basis(k))

# Basis checking
check_basis(::Union{Ket{B1}, WeightedKet{B1}, SumKet{B1}}, ::Union{Ket{B2}, WeightedKet{B2}, SumKet{B2}}) where {B1,B2} =
    B1 == B2 || throw(DimensionMismatch("Kets are in different bases"))
check_basis(::Union{Bra{B1}, WeightedBra{B1}, SumBra{B1}}, ::Union{Bra{B2}, WeightedBra{B2}, SumBra{B2}}) where {B1,B2} =
    B1 == B2 || throw(DimensionMismatch("Bras are in different bases"))

check_space(k1::AbstractKet, k2::AbstractKet) = 
    space(k1) == space(k2) || throw(DimensionMismatch("Kets are in different spaces"))
check_space(b1::AbstractBra, b2::AbstractBra) = 
    space(b1) == space(b2) || throw(DimensionMismatch("Bras are in different spaces"))

# Equality
Base.:(==)(k1::Ket{B1}, k2::Ket{B2}) where {B1, B2} = B1 == B2 && k1.index == k2.index
Base.:(==)(b1::Bra{B1}, b2::Bra{B2}) where {B1, B2} = B1 == B2 && b1.index == b2.index
Base.:(==)(pk1::ProductKet{A1,A2}, pk2::ProductKet{B1,B2}) where {A1,A2,B1,B2} = 
    A1 == B1 && A2 == B2 && pk1.ket1 == pk2.ket1 && pk1.ket2 == pk2.ket2
Base.:(==)(pb1::ProductBra{A1,A2}, pb2::ProductBra{B1,B2}) where {A1,A2,B1,B2} = 
    A1 == B1 && A2 == B2 && pb1.bra1 == pb2.bra1 && pb1.bra2 == pb2.bra2

# Helper to format index (single or multi)
function _format_index(idx)
    if isnothing(idx)
        "ψ"
    elseif idx isa Tuple
        join(string.(idx), ",")
    else
        string(idx)
    end
end

# Show methods
function Base.show(io::IO, k::Ket)
    name = _format_index(k.index)
    print(io, "|", name, "⟩")
end

function Base.show(io::IO, b::Bra)
    name = _format_index(b.index)
    print(io, "⟨", name, "|")
end

function Base.show(io::IO, pk::ProductKet)
    print(io, pk.ket1, "⊗", pk.ket2)
end

function Base.show(io::IO, pb::ProductBra)
    print(io, pb.bra1, "⊗", pb.bra2)
end

function Base.show(io::IO, wk::WeightedKet)
    print(io, wk.weight, "·")
    if wk.ket isa Ket
        print(io, wk.ket)
    else  # ProductKet
        print(io, wk.ket)
    end
end

function Base.show(io::IO, wb::WeightedBra)
    print(io, wb.weight, "·")
    if wb.bra isa Bra
        print(io, wb.bra)
    else  # ProductBra
        print(io, wb.bra)
    end
end

function Base.show(io::IO, sk::SumKet)
    if !isnothing(sk.display_name)
        print(io, "|", sk.display_name, "⟩")
    else
        for (i, (k, w)) in enumerate(zip(sk.kets, sk.weights))
            i > 1 && print(io, w >= 0 ? " + " : " - ")
            i == 1 && w < 0 && print(io, "-")
            abs(w) != 1 && print(io, abs(w), "·")
            print(io, k)
        end
    end
end

function Base.show(io::IO, sb::SumBra)
    if !isnothing(sb.display_name)
        print(io, "⟨", sb.display_name, "|")
    else
        for (i, (b, w)) in enumerate(zip(sb.bras, sb.weights))
            i > 1 && print(io, w >= 0 ? " + " : " - ")
            i == 1 && w < 0 && print(io, "-")
            abs(w) != 1 && print(io, abs(w), "·")
            print(io, b)
        end
    end
end

# Convenience constructors (Fock space)
@doc """
    FockKet(space::HilbertSpace, n::Int)

Create a Fock state |n⟩ in the given infinite-dimensional Hilbert space.

See also: [`Ket`](@ref), [`FockBra`](@ref)
""" FockKet
function FockKet(space::HilbertSpace{T,dim}, n::Int) where {T,dim}
    dim isa Tuple{Nothing} || throw(ArgumentError("Not a valid Fock space (dimension is limited)"))
    Ket(space, n)
end

@doc """
    FockBra(space::HilbertSpace, n::Int)

Create a Fock bra ⟨n| in the given infinite-dimensional Hilbert space.

See also: [`Bra`](@ref), [`FockKet`](@ref)
""" FockBra
function FockBra(space::HilbertSpace{T,dim}, n::Int) where {T,dim}
    dim isa Tuple{Nothing} || throw(ArgumentError("Not a valid Fock space (dimension is limited)"))
    Bra(space, n)
end
