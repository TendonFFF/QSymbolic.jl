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
    ProductKet{Bs<:Tuple}(kets::Vector{Ket})
    ProductKet(ket1, ket2, ...)

Tensor product of kets: |ψ₁⟩⊗|ψ₂⟩⊗...⊗|ψₙ⟩.
Lives in a `CompositeBasis`. This is a container struct (like SumKet),
not a basic element - only `Ket` is a basic element.

**Note**: ProductKets are **order-independent** (bosonic/symmetric).
The kets are canonically ordered by basis, so `k1 ⊗ k2 == k2 ⊗ k1` when
they have different bases. This is the default behavior; fermionic
(anti-symmetric) tensor products will be added in a future update.

The tensor product operator `⊗` creates ProductKets automatically.

# Examples
```jldoctest
julia> H1, H2 = HilbertSpace(:A, 2), HilbertSpace(:B, 2);

julia> ψ, ϕ = Ket(H1, :ψ), Ket(H2, :ϕ);

julia> ψ ⊗ ϕ
|ψ⟩⊗|ϕ⟩

julia> ϕ ⊗ ψ == ψ ⊗ ϕ  # Order-independent (bosonic)
true

julia> H3 = HilbertSpace(:C, 2);

julia> χ = Ket(H3, :χ);

julia> ψ ⊗ ϕ ⊗ χ  # Three kets
|ψ⟩⊗|ϕ⟩⊗|χ⟩
```

See also: [`Ket`](@ref), [`⊗`](@ref)
""" ProductKet
struct ProductKet{Bs<:Tuple} <: AbstractKet{CompositeBasis{Bs}}
    kets::Vector{Ket}  # Vector of basic kets (canonically ordered)
    
    function ProductKet(kets::Vector{Ket})
        if length(kets) < 2
            throw(ArgumentError("ProductKet requires at least 2 kets"))
        end
        # Canonicalize order: sort kets by basis type for order-independence
        sorted_kets = _canonicalize_product_kets(kets)
        # Extract basis types from each ket
        basis_types = Tuple{[basis(typeof(k)) for k in sorted_kets]...}
        new{basis_types}(sorted_kets)
    end
    
    # Convenience constructor from varargs
    function ProductKet(ket1::Ket, ket2::Ket, kets::Ket...)
        ProductKet([ket1, ket2, kets...])
    end
end

# Canonicalize ket order for bosonic (order-independent) behavior
function _canonicalize_product_kets(kets::Vector{Ket})
    # Sort by basis type, then by index for deterministic ordering
    return sort(kets, by = k -> (string(basis(typeof(k))), string(k.index)))
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
    function WeightedKet(pk::ProductKet{Bs}, weight::T) where {Bs,T}
        CB = CompositeBasis{Bs}
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
    function SumKet(ket::Union{Ket{B}, ProductKet{Bs}}; name::Union{Symbol, Nothing}=nothing) where {B, Bs}
        if ket isa Ket
            new{B, Int}(name, [ket], [1])
        else  # ProductKet
            CB = CompositeBasis{Bs}
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

# ==================== UTILITY FUNCTIONS ====================

# basis() and space() functions
basis(::Ket{B}) where B = B
basis(::Bra{B}) where B = B
basis(::Type{Ket{B}}) where B = B
basis(::Type{Bra{B}}) where B = B
basis(::ProductKet{Bs}) where {Bs} = CompositeBasis{Bs}
basis(::ProductBra{Bs}) where {Bs} = CompositeBasis{Bs}
basis(::Type{ProductKet{Bs}}) where {Bs} = CompositeBasis{Bs}
basis(::Type{ProductBra{Bs}}) where {Bs} = CompositeBasis{Bs}
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

# ProductKet/ProductBra equality: order-independent (bosonic)
# Since kets/bras are stored in canonical order, simple comparison works
Base.:(==)(pk1::ProductKet{Bs1}, pk2::ProductKet{Bs2}) where {Bs1,Bs2} = 
    Bs1 == Bs2 && all(pk1.kets[i] == pk2.kets[i] for i in 1:length(pk1.kets))
Base.:(==)(pb1::ProductBra{Bs1}, pb2::ProductBra{Bs2}) where {Bs1,Bs2} = 
    Bs1 == Bs2 && all(pb1.bras[i] == pb2.bras[i] for i in 1:length(pb1.bras))

# Hash functions for proper collection behavior
Base.hash(k::Ket, h::UInt) = hash((basis(typeof(k)), k.index), h)
Base.hash(b::Bra, h::UInt) = hash((basis(typeof(b)), b.index), h)
Base.hash(pk::ProductKet, h::UInt) = hash((basis(typeof(pk)), pk.kets), h)
Base.hash(pb::ProductBra, h::UInt) = hash((basis(typeof(pb)), pb.bras), h)

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
    for (i, k) in enumerate(pk.kets)
        i > 1 && print(io, "⊗")
        print(io, k)
    end
end

function Base.show(io::IO, pb::ProductBra)
    for (i, b) in enumerate(pb.bras)
        i > 1 && print(io, "⊗")
        print(io, b)
    end
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
