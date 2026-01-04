# Basis types for Hilbert spaces

@doc """
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
""" Basis
struct Basis{S<:AbstractSpace, name} <: AbstractBasis{S}
    function Basis(space::S, name::Symbol) where S<:AbstractSpace
        new{S, name}()
    end
    function Basis{S,name}() where {S<:AbstractSpace, name}
        new{S, name}()
    end
end

@doc """
    space(basis)
    space(ket_or_bra)
    space(operator)

Get the underlying Hilbert space of a basis, ket, bra, or operator.
""" space
space(::Type{Basis{S,name}}) where {S,name} = S
space(::Basis{S,name}) where {S,name} = S

@doc """
    basisname(basis)

Get the symbolic name of a basis.
""" basisname
basisname(::Type{Basis{S,name}}) where {S,name} = name
basisname(::Basis{S,name}) where {S,name} = name

@eval Base.$(:(==))(::Basis{S1,n1}, ::Basis{S2,n2}) where {S1,S2,n1,n2} = S1 == S2 && n1 == n2

@doc """
    CompositeBasis{Bs<:Tuple}

Tensor product of bases. Created via `basis1 ⊗ basis2` or `basis1 ⊗ basis2 ⊗ basis3 ⊗ ...`.

Supports arbitrary number of component bases using tuple type parameters.

Transforms factorize automatically: if transforms exist for each component basis,
tensor product transforms are automatically derived without explicit definition.

# Examples
```jldoctest
julia> H1, H2 = HilbertSpace(:A, 2), HilbertSpace(:B, 2);

julia> Za, Zb = Basis(H1, :z), Basis(H2, :z);

julia> Za ⊗ Zb
Basis{z}⊗Basis{z}

julia> H3 = HilbertSpace(:C, 2);

julia> Zc = Basis(H3, :z);

julia> Za ⊗ Zb ⊗ Zc  # Three bases
Basis{z}⊗Basis{z}⊗Basis{z}
```

See also: [`Basis`](@ref), [`define_transform!`](@ref)
""" CompositeBasis
struct CompositeBasis{Bs<:Tuple} <: AbstractBasis{CompositeSpace}
    function CompositeBasis{Bs}() where {Bs<:Tuple}
        # Validate all elements are AbstractBasis types
        for B in Bs.parameters
            B <: AbstractBasis || throw(ArgumentError("All elements must be AbstractBasis types"))
        end
        new{Bs}()
    end
end

# Constructor for tuple of basis types
CompositeBasis(bases::Type{<:AbstractBasis}...) = CompositeBasis{Tuple{bases...}}()

@doc """
    ⊗(basis1::AbstractBasis, basis2::AbstractBasis)
    ⊗(cb::CompositeBasis, basis::AbstractBasis)
    ⊗(basis::AbstractBasis, cb::CompositeBasis)

Tensor product of bases. Supports chaining for multiple bases.
"""
function ⊗(::B1, ::B2) where {B1<:AbstractBasis, B2<:AbstractBasis}
    # Create composite basis with two elements
    CompositeBasis{Tuple{B1, B2}}()
end

# Allow chaining: CompositeBasis ⊗ Basis
function ⊗(::CompositeBasis{Bs}, ::B) where {Bs, B<:AbstractBasis}
    # Add new basis to the tuple
    CompositeBasis{Tuple{Bs.parameters..., B}}()
end

# Allow chaining: Basis ⊗ CompositeBasis
function ⊗(::B, ::CompositeBasis{Bs}) where {B<:AbstractBasis, Bs}
    # Prepend basis to the tuple
    CompositeBasis{Tuple{B, Bs.parameters...}}()
end

# CompositeBasis ⊗ CompositeBasis
function ⊗(::CompositeBasis{Bs1}, ::CompositeBasis{Bs2}) where {Bs1, Bs2}
    # Concatenate tuples
    CompositeBasis{Tuple{Bs1.parameters..., Bs2.parameters...}}()
end

@doc """
    basis(ket_or_bra)
    basis(operator)

Get the basis type of a ket, bra, or operator.
""" basis

@doc """
    basis1(cb::CompositeBasis)

Get the first component basis of a composite basis.
""" basis1
basis1(::Type{CompositeBasis{Bs}}) where {Bs} = Bs.parameters[1]
basis1(::CompositeBasis{Bs}) where {Bs} = Bs.parameters[1]

@doc """
    basis2(cb::CompositeBasis)

Get the second component basis of a composite basis.
""" basis2
basis2(::Type{CompositeBasis{Bs}}) where {Bs} = Bs.parameters[2]
basis2(::CompositeBasis{Bs}) where {Bs} = Bs.parameters[2]

@doc """
    bases(cb::CompositeBasis)

Get all component bases of a composite basis as a tuple.
""" bases
bases(::Type{CompositeBasis{Bs}}) where {Bs} = Bs.parameters
bases(::CompositeBasis{Bs}) where {Bs} = Bs.parameters

# Space of a composite basis is the composite of its component spaces
function _composite_space_type(basis_types::Tuple)
    # Extract space for each basis type
    spaces = [space(B) for B in basis_types]
    
    # Flatten all space parameters into tuples
    # Collect all names and dims from each space
    names_vec = vcat([collect(S.parameters[1]) for S in spaces]...)
    dims_vec = vcat([collect(S.parameters[2]) for S in spaces]...)
    all_names = Tuple(names_vec)
    all_dims = Tuple(dims_vec)
    
    return CompositeSpace{all_names, all_dims}
end

space(::Type{CompositeBasis{Bs}}) where {Bs} = _composite_space_type(Bs.parameters)
space(::CompositeBasis{Bs}) where {Bs} = _composite_space_type(Bs.parameters)

@eval Base.$(:(==))(::CompositeBasis{Bs1}, ::CompositeBasis{Bs2}) where {Bs1,Bs2} = Bs1 == Bs2

@doc """
    DefaultBasis{S}

Implicit basis when no explicit basis is specified. Provides backward compatibility
for code that creates kets directly from spaces without specifying a basis.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2);

julia> ψ = BasisKet(H, :ψ)  # uses DefaultBasis internally
|ψ⟩
```
""" DefaultBasis
struct DefaultBasis{S<:AbstractSpace} <: AbstractBasis{S} end

space(::Type{DefaultBasis{S}}) where S = S
space(::DefaultBasis{S}) where S = S

# Show methods
function Base.show(io::IO, b::Basis{S,name}) where {S,name}
    print(io, "Basis{", name, "}")
end

function Base.show(io::IO, ::CompositeBasis{Bs}) where {Bs}
    # Show each basis component separated by ⊗
    for (i, B) in enumerate(Bs.parameters)
        i > 1 && print(io, "⊗")
        show(io, B())
    end
end

# Iteration protocol for HilbertSpace to enable destructuring
# Allows: H, Hb = HilbertSpace(:H, 2) to get both space and default basis
# Defined here after Basis type is available
@doc """
    iterate(space::HilbertSpace)
    iterate(space::HilbertSpace, state)

Iteration protocol for `HilbertSpace` that enables convenient destructuring.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2)  # Just the space
ℋ(H, dim=2)

julia> H, Hb = HilbertSpace(:H, 2)  # Space and default basis
(ℋ(H, dim=2), Basis{default})
```

The first iteration returns the space itself, and the second returns a default basis.
This allows for convenient syntax where you can get both the space and a default basis in one line.
""" Base.iterate
function Base.iterate(space::HilbertSpace)
    # First iteration: return the space itself, with state=1
    return (space, 1)
end

function Base.iterate(space::HilbertSpace, state::Int)
    if state == 1
        # Second iteration: return a default basis, with state=2
        return (Basis(space, :default), 2)
    else
        # No more iterations
        return nothing
    end
end

# Declare HilbertSpace as iterable with length 2
Base.length(::HilbertSpace) = 2
