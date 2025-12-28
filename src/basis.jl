# Basis types for Hilbert spaces

"""
    Basis(space, name::Symbol)

A named orthonormal basis for a Hilbert space. Kets in same basis
are orthonormal; cross-basis inner products require a transform.

# Examples
```jldoctest
julia> H = HilbertSpace(:spin, 2);

julia> Zb = Basis(H, :z)  # spin-z basis
Basis{z}

julia> Xb = Basis(H, :x)  # spin-x basis
Basis{x}
```

See also: [`define_transform!`](@ref), [`DefaultBasis`](@ref), [`CompositeBasis`](@ref)
"""
struct Basis{S<:AbstractSpace, name} <: AbstractBasis{S}
    function Basis(space::S, name::Symbol) where S<:AbstractSpace
        new{S, name}()
    end
    function Basis{S,name}() where {S<:AbstractSpace, name}
        new{S, name}()
    end
end

"""
    space(basis)

Get the underlying Hilbert space of a basis.
"""
space(::Type{Basis{S,name}}) where {S,name} = S
space(::Basis{S,name}) where {S,name} = S

"""
    basisname(basis)

Get the symbolic name of a basis.
"""
basisname(::Type{Basis{S,name}}) where {S,name} = name
basisname(::Basis{S,name}) where {S,name} = name

@eval Base.$(:(==))(::Basis{S1,n1}, ::Basis{S2,n2}) where {S1,S2,n1,n2} = S1 == S2 && n1 == n2

"""
    CompositeBasis{B1,B2}

Tensor product of two bases. Created via `basis1 ⊗ basis2`.

Transforms factorize automatically: if transforms exist for `B1→B1'` and `B2→B2'`,
then `B1⊗B2 → B1'⊗B2'` is automatically derived without explicit definition.

# Examples
```jldoctest
julia> H1, H2 = HilbertSpace(:A, 2), HilbertSpace(:B, 2);

julia> Za, Zb = Basis(H1, :z), Basis(H2, :z);

julia> Za ⊗ Zb
Basis{z}⊗Basis{z}
```

See also: [`Basis`](@ref), [`define_transform!`](@ref)
"""
struct CompositeBasis{B1<:AbstractBasis, B2<:AbstractBasis} <: AbstractBasis{CompositeSpace}
    function CompositeBasis{B1,B2}() where {B1<:AbstractBasis, B2<:AbstractBasis}
        new{B1,B2}()
    end
end

"""
    ⊗(basis1::AbstractBasis, basis2::AbstractBasis)

Tensor product of two bases.
"""
⊗(::B1, ::B2) where {B1<:AbstractBasis, B2<:AbstractBasis} = CompositeBasis{B1,B2}()

"""
    basis1(cb::CompositeBasis)

Get the first component basis of a composite basis.
"""
basis1(::Type{CompositeBasis{B1,B2}}) where {B1,B2} = B1
basis1(::CompositeBasis{B1,B2}) where {B1,B2} = B1

"""
    basis2(cb::CompositeBasis)

Get the second component basis of a composite basis.
"""
basis2(::Type{CompositeBasis{B1,B2}}) where {B1,B2} = B2
basis2(::CompositeBasis{B1,B2}) where {B1,B2} = B2

space(::Type{CompositeBasis{B1,B2}}) where {B1,B2} = CompositeSpace
space(::CompositeBasis{B1,B2}) where {B1,B2} = CompositeSpace

@eval Base.$(:(==))(::CompositeBasis{A1,A2}, ::CompositeBasis{B1,B2}) where {A1,A2,B1,B2} = A1 == B1 && A2 == B2

"""
    DefaultBasis{S}

Implicit basis when no explicit basis is specified. Provides backward compatibility
for code that creates kets directly from spaces without specifying a basis.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ = BasisKet(H, :ψ)  # uses DefaultBasis internally
|ψ⟩
```
"""
struct DefaultBasis{S<:AbstractSpace} <: AbstractBasis{S} end

space(::Type{DefaultBasis{S}}) where S = S
space(::DefaultBasis{S}) where S = S
