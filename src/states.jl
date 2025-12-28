# Ket and Bra state vectors

"""
    BasisKet(basis, index)
    BasisKet(space, index)  # uses DefaultBasis

A basis ket |index⟩ in the specified basis.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> Zb = Basis(H, :z);

julia> up = BasisKet(Zb, :↑)   # |↑⟩ in z-basis
|↑⟩

julia> ψ = BasisKet(H, :ψ)     # |ψ⟩ in default basis
|ψ⟩
```

See also: [`BasisBra`](@ref), [`weightedKet`](@ref), [`sumKet`](@ref)
"""
struct BasisKet{B<:AbstractBasis} <: AbstractKet{B}
    index::Union{Symbol, Nothing}

    function BasisKet(::B, idx::Union{Symbol, Nothing}=nothing) where B<:AbstractBasis
        new{B}(idx)
    end
    function BasisKet(::B, idx::Int) where B<:AbstractBasis
        new{B}(Symbol(idx))
    end
    function BasisKet{B}(idx::Union{Symbol, Int, Nothing}=nothing) where B<:AbstractBasis
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function BasisKet(space::AbstractSpace, idx::Union{Symbol, Int, Nothing}=nothing)
        B = DefaultBasis{typeof(space)}
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
end

"""
    weightedKet{B,T}

Ket multiplied by scalar weight. Created via `weight * ket`.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ = BasisKet(H, :ψ);

julia> 2 * ψ
2·|ψ⟩
```
"""
struct weightedKet{B<:AbstractBasis, T} <: AbstractKet{B}
    Ket::BasisKet{B}
    weight::T
end

"""
    sumKet{B,T}

Linear combination of basis kets. Created via ket addition.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ, ϕ = BasisKet(H, :ψ), BasisKet(H, :ϕ);

julia> ψ + ϕ
|ψ⟩ + |ϕ⟩
```
"""
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

"""
    BasisBra(basis, index)
    BasisBra(space, index)  # uses DefaultBasis

A basis bra ⟨index| in the specified basis. Also created via `ket'`.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ = BasisKet(H, :ψ);

julia> ψ'
⟨ψ|
```
"""
struct BasisBra{B<:AbstractBasis} <: AbstractBra{B}
    index::Union{Symbol, Nothing}

    function BasisBra(::B, idx::Union{Symbol, Nothing}=nothing) where B<:AbstractBasis
        new{B}(idx)
    end
    function BasisBra(::B, idx::Int) where B<:AbstractBasis
        new{B}(Symbol(idx))
    end
    function BasisBra{B}(idx::Union{Symbol, Int, Nothing}=nothing) where B<:AbstractBasis
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
    function BasisBra(space::AbstractSpace, idx::Union{Symbol, Int, Nothing}=nothing)
        B = DefaultBasis{typeof(space)}
        new{B}(isnothing(idx) ? nothing : Symbol(idx))
    end
end

"""
    weightedBra{B,T}

Bra multiplied by scalar weight. Created via `weight * bra`.
"""
struct weightedBra{B<:AbstractBasis, T} <: AbstractBra{B}
    Bra::BasisBra{B}
    weight::T
end

"""
    sumBra{B,T}

Linear combination of basis bras. Created via bra addition.
"""
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

"""
    ProductKet{B1,B2}(ket1, ket2)

Tensor product of two kets: |ψ⟩⊗|ϕ⟩. Lives in `CompositeBasis{B1,B2}`.

# Examples
```jldoctest
julia> H1, H2 = HilbertSpace(:A, 2), HilbertSpace(:B, 2);

julia> ψ, ϕ = BasisKet(H1, :ψ), BasisKet(H2, :ϕ);

julia> ψ ⊗ ϕ
|ψ⟩⊗|ϕ⟩
```
"""
struct ProductKet{B1<:AbstractBasis, B2<:AbstractBasis} <: AbstractKet{CompositeBasis{B1,B2}}
    ket1::BasisKet{B1}
    ket2::BasisKet{B2}
end

"""
    ProductBra{B1,B2}(bra1, bra2)

Tensor product of two bras: ⟨ψ|⊗⟨ϕ|. Lives in `CompositeBasis{B1,B2}`.
"""
struct ProductBra{B1<:AbstractBasis, B2<:AbstractBasis} <: AbstractBra{CompositeBasis{B1,B2}}
    bra1::BasisBra{B1}
    bra2::BasisBra{B2}
end

"""
    SumProductKet{B1,B2,T}

Linear combination of product kets. Used for entangled states and 
transformed composite states.

# Examples
```julia
# Bell state |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
bell = (zero_zero + one_one) / √2
```
"""
struct SumProductKet{B1<:AbstractBasis, B2<:AbstractBasis, T} <: AbstractKet{CompositeBasis{B1,B2}}
    display_name::Union{Symbol,Nothing}
    kets::Vector{ProductKet{B1,B2}}
    weights::Vector{T}
    
    function SumProductKet(kets::AbstractVector{ProductKet{B1,B2}}, weights::AbstractVector{T}; name::Union{Symbol,Nothing}=nothing) where {B1,B2,T}
        new{B1,B2,T}(name, kets, weights)
    end
end

"""
    SumProductBra{B1,B2,T}

Linear combination of product bras.
"""
struct SumProductBra{B1<:AbstractBasis, B2<:AbstractBasis, T} <: AbstractBra{CompositeBasis{B1,B2}}
    display_name::Union{Symbol,Nothing}
    bras::Vector{ProductBra{B1,B2}}
    weights::Vector{T}
    
    function SumProductBra(bras::AbstractVector{ProductBra{B1,B2}}, weights::AbstractVector{T}; name::Union{Symbol,Nothing}=nothing) where {B1,B2,T}
        new{B1,B2,T}(name, bras, weights)
    end
end

"""
    InnerProduct{B1,B2}(bra, ket)

Symbolic inner product ⟨bra|ket⟩ for cross-basis when no transform is defined.
Returned when computing inner products between different bases without a registered transform.

See also: [`define_transform!`](@ref)
"""
struct InnerProduct{B1<:AbstractBasis, B2<:AbstractBasis}
    bra::BasisBra{B1}
    ket::BasisKet{B2}
end

# Utility functions

"""
    basis(ket_or_bra)

Get the basis type of a ket or bra.
"""
basis(::AbstractKet{B}) where B = B
basis(::AbstractBra{B}) where B = B
basis(::Type{<:AbstractKet{B}}) where B = B
basis(::Type{<:AbstractBra{B}}) where B = B

"""
    space(ket_or_bra)

Get the underlying Hilbert space of a ket or bra.
"""
space(::AbstractKet{B}) where B = space(B)
space(::AbstractBra{B}) where B = space(B)
