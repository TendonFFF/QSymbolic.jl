# Ket and Bra state vectors

# Single index type - supports symbolic indices
const SingleIndex = Union{Symbol, Nothing, AbstractSymbolic}

# Multi-index type - tuple of single indices (for composite systems)
const MultiIndex = Tuple{Vararg{SingleIndex}}

# Combined index type - single or multi-index
const KetIndex = Union{SingleIndex, MultiIndex}

@doc """
    BasisKet(basis, index)
    BasisKet(space, index)  # uses DefaultBasis

A basis ket |index⟩ in the specified basis. 

Index can be:
- A `Symbol` (e.g., `:↑`, `:ψ`)
- An `Int` (converted to Symbol)
- A symbolic expression (`Sym`, `SymExpr`)
- A tuple for multi-index (e.g., `(:n, :m)` for composite bases)

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> Zb = Basis(H, :z);

julia> up = BasisKet(Zb, :↑)   # |↑⟩ in z-basis
|↑⟩

julia> ψ = BasisKet(H, :ψ)     # |ψ⟩ in default basis
|ψ⟩
```

Symbolic indices are also supported:
```julia
n = Sym(:n)
ket_n = BasisKet(Fock_basis, n)  # |n⟩ with symbolic index
```

Multi-indices for composite bases:
```julia
composite_basis = B1 ⊗ B2
ket = BasisKet(composite_basis, (Sym(:n), Sym(:m)))  # |n,m⟩
```

See also: [`BasisBra`](@ref), [`weightedKet`](@ref), [`sumKet`](@ref)
""" BasisKet
struct BasisKet{B<:AbstractBasis} <: AbstractKet{B}
    index::KetIndex

    function BasisKet(::B, idx::KetIndex=nothing) where B<:AbstractBasis
        new{B}(idx)
    end
    function BasisKet(::B, idx::Int) where B<:AbstractBasis
        new{B}(Symbol(idx))
    end
    # Multi-index from tuple of integers
    function BasisKet(::B, idx::Tuple{Vararg{Int}}) where B<:AbstractBasis
        new{B}(Symbol.(idx))
    end
    function BasisKet{B}(idx::Union{Symbol, Int, Nothing}=nothing) where B<:AbstractBasis
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function BasisKet{B}(idx::AbstractSymbolic) where B<:AbstractBasis
        new{B}(idx)
    end
    function BasisKet{B}(idx::MultiIndex) where B<:AbstractBasis
        new{B}(idx)
    end
    function BasisKet(space::AbstractSpace, idx::Union{Symbol, Int, Nothing}=nothing)
        B = DefaultBasis{typeof(space)}
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function BasisKet(space::AbstractSpace, idx::AbstractSymbolic)
        B = DefaultBasis{typeof(space)}
        new{B}(idx)
    end
    function BasisKet(space::AbstractSpace, idx::MultiIndex)
        B = DefaultBasis{typeof(space)}
        new{B}(idx)
    end
    # Multi-index from tuple of integers (space version)
    function BasisKet(space::AbstractSpace, idx::Tuple{Vararg{Int}})
        B = DefaultBasis{typeof(space)}
        new{B}(Symbol.(idx))
    end
end

@doc """
    weightedKet{B,T}

Ket multiplied by scalar weight. Created via `weight * ket`.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ = BasisKet(H, :ψ);

julia> 2 * ψ
2·|ψ⟩
```
""" weightedKet
struct weightedKet{B<:AbstractBasis, T} <: AbstractKet{B}
    Ket::BasisKet{B}
    weight::T
end

@doc """
    sumKet{B,T}

Linear combination of basis kets. Created via ket addition.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ, ϕ = BasisKet(H, :ψ), BasisKet(H, :ϕ);

julia> ψ + ϕ
|ψ⟩ + |ϕ⟩
```
""" sumKet
struct sumKet{B<:AbstractBasis, T} <: AbstractKet{B}
    display_name::Union{Symbol,Nothing}
    kets::Vector{BasisKet{B}}
    weights::Vector{T}

    function sumKet(kets::AbstractVector{BasisKet{B}}, weights::AbstractVector{T}; name::Union{Symbol, Nothing}=nothing) where {B, T}
        new{B, T}(name, kets, weights)
    end
    function sumKet(Ket::BasisKet{B}; name::Union{Symbol, Nothing}=nothing) where B
        new{B, Int}(name, [Ket], [1])
    end
    function sumKet(ket::weightedKet{B, T}; name::Union{Symbol, Nothing}=nothing) where {B, T}
        new{B, T}(name, [ket.Ket], [ket.weight])
    end
end

@doc """
    BasisBra(basis, index)
    BasisBra(space, index)  # uses DefaultBasis

A basis bra ⟨index| in the specified basis. Also created via `ket'`.

Index can be:
- A `Symbol` (e.g., `:↑`, `:ψ`)
- An `Int` (converted to Symbol)
- A symbolic expression (`Sym`, `SymExpr`)
- A tuple for multi-index (e.g., `(:n, :m)` for composite bases)

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ = BasisKet(H, :ψ);

julia> ψ'
⟨ψ|
```
""" BasisBra
struct BasisBra{B<:AbstractBasis} <: AbstractBra{B}
    index::KetIndex

    function BasisBra(::B, idx::KetIndex=nothing) where B<:AbstractBasis
        new{B}(idx)
    end
    function BasisBra(::B, idx::Int) where B<:AbstractBasis
        new{B}(Symbol(idx))
    end
    # Multi-index from tuple of integers
    function BasisBra(::B, idx::Tuple{Vararg{Int}}) where B<:AbstractBasis
        new{B}(Symbol.(idx))
    end
    function BasisBra{B}(idx::Union{Symbol, Int, Nothing}=nothing) where B<:AbstractBasis
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function BasisBra{B}(idx::AbstractSymbolic) where B<:AbstractBasis
        new{B}(idx)
    end
    function BasisBra{B}(idx::MultiIndex) where B<:AbstractBasis
        new{B}(idx)
    end
    function BasisBra(space::AbstractSpace, idx::Union{Symbol, Int, Nothing}=nothing)
        B = DefaultBasis{typeof(space)}
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function BasisBra(space::AbstractSpace, idx::AbstractSymbolic)
        B = DefaultBasis{typeof(space)}
        new{B}(idx)
    end
    function BasisBra(space::AbstractSpace, idx::MultiIndex)
        B = DefaultBasis{typeof(space)}
        new{B}(idx)
    end
    # Multi-index from tuple of integers (space version)
    function BasisBra(space::AbstractSpace, idx::Tuple{Vararg{Int}})
        B = DefaultBasis{typeof(space)}
        new{B}(Symbol.(idx))
    end
end

@doc """
    weightedBra{B,T}

Bra multiplied by scalar weight. Created via `weight * bra`.
""" weightedBra
struct weightedBra{B<:AbstractBasis, T} <: AbstractBra{B}
    Bra::BasisBra{B}
    weight::T
end

@doc """
    sumBra{B,T}

Linear combination of basis bras. Created via bra addition.
""" sumBra
struct sumBra{B<:AbstractBasis, T} <: AbstractBra{B}
    display_name::Union{Symbol,Nothing}
    bras::Vector{BasisBra{B}}
    weights::Vector{T}

    function sumBra(bras::AbstractVector{BasisBra{B}}, weights::AbstractVector{T}; name::Union{Symbol, Nothing}=nothing) where {B, T}
        new{B, T}(name, bras, weights)
    end
    function sumBra(Bra::BasisBra{B}; name::Union{Symbol, Nothing}=nothing) where B
        new{B, Int}(name, [Bra], [1])
    end
    function sumBra(bra::weightedBra{B, T}; name::Union{Symbol, Nothing}=nothing) where {B, T}
        new{B, T}(name, [bra.Bra], [bra.weight])
    end
end

# Product states for composite systems

@doc """
    ProductKet{B1,B2}(ket1, ket2)

Tensor product of two kets: |ψ⟩⊗|ϕ⟩. Lives in `CompositeBasis{B1,B2}`.

# Examples
```jldoctest
julia> H1, H2 = HilbertSpace(:A, 2), HilbertSpace(:B, 2);

julia> ψ, ϕ = BasisKet(H1, :ψ), BasisKet(H2, :ϕ);

julia> ψ ⊗ ϕ
|ψ⟩⊗|ϕ⟩
```
""" ProductKet
struct ProductKet{B1<:AbstractBasis, B2<:AbstractBasis} <: AbstractKet{CompositeBasis{B1,B2}}
    ket1::BasisKet{B1}
    ket2::BasisKet{B2}
end

@doc """
    ProductBra{B1,B2}(bra1, bra2)

Tensor product of two bras: ⟨ψ|⊗⟨ϕ|. Lives in `CompositeBasis{B1,B2}`.
""" ProductBra
struct ProductBra{B1<:AbstractBasis, B2<:AbstractBasis} <: AbstractBra{CompositeBasis{B1,B2}}
    bra1::BasisBra{B1}
    bra2::BasisBra{B2}
end

@doc """
    SumProductKet{B1,B2,T}

Linear combination of product kets. Used for entangled states and 
transformed composite states.

# Examples
```julia
# Bell state |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
bell = (zero_zero + one_one) / √2
```
""" SumProductKet
struct SumProductKet{B1<:AbstractBasis, B2<:AbstractBasis, T} <: AbstractKet{CompositeBasis{B1,B2}}
    display_name::Union{Symbol,Nothing}
    kets::Vector{ProductKet{B1,B2}}
    weights::Vector{T}
    
    function SumProductKet(kets::AbstractVector{ProductKet{B1,B2}}, weights::AbstractVector{T}; name::Union{Symbol,Nothing}=nothing) where {B1,B2,T}
        new{B1,B2,T}(name, kets, weights)
    end
end

@doc """
    SumProductBra{B1,B2,T}

Linear combination of product bras.
""" SumProductBra
struct SumProductBra{B1<:AbstractBasis, B2<:AbstractBasis, T} <: AbstractBra{CompositeBasis{B1,B2}}
    display_name::Union{Symbol,Nothing}
    bras::Vector{ProductBra{B1,B2}}
    weights::Vector{T}
    
    function SumProductBra(bras::AbstractVector{ProductBra{B1,B2}}, weights::AbstractVector{T}; name::Union{Symbol,Nothing}=nothing) where {B1,B2,T}
        new{B1,B2,T}(name, bras, weights)
    end
end

@doc """
    InnerProduct{B1,B2}(bra, ket)

Symbolic inner product ⟨bra|ket⟩ for cross-basis when no transform is defined.
Returned when computing inner products between different bases without a registered transform.

See also: [`define_transform!`](@ref)
""" InnerProduct
struct InnerProduct{B1<:AbstractBasis, B2<:AbstractBasis}
    bra::BasisBra{B1}
    ket::BasisKet{B2}
end

# Utility functions

@doc """
    basis(ket_or_bra)

Get the basis type of a ket or bra.
""" basis
basis(::AbstractKet{B}) where B = B
basis(::AbstractBra{B}) where B = B
basis(::Type{<:AbstractKet{B}}) where B = B
basis(::Type{<:AbstractBra{B}}) where B = B

@doc """
    space(ket_or_bra)

Get the underlying Hilbert space of a ket or bra.
""" space
space(::AbstractKet{B}) where B = space(B)
space(::AbstractBra{B}) where B = space(B)
