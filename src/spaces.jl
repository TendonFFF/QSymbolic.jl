# Hilbert spaces and composite spaces

"""
    HilbertSpace(name::Symbol, dim=nothing)

A named Hilbert space with optional finite dimension.

# Examples
```jldoctest
julia> H = HilbertSpace(:H, 2)  # qubit
ℋ(H, dim=2)

julia> F = HilbertSpace(:F)     # infinite-dim
ℋ(F)
```
"""
struct HilbertSpace{name, dim} <: AbstractSpace{name, dim}
    function HilbertSpace(A::Symbol, dim::Union{Int,Nothing}=nothing)
        new{(A,), (dim,)}()
    end
end

"""
    CompositeSpace(space1, space2)

Tensor product of two spaces. Usually created via `space1 ⊗ space2`.

# Examples
```jldoctest
julia> H1 = HilbertSpace(:A, 2);

julia> H2 = HilbertSpace(:B, 3);

julia> H1 ⊗ H2
ℋ(A) ⊗ ℋ(B)
```
"""
struct CompositeSpace{name, dim} <: AbstractSpace{name, dim}
    function CompositeSpace(space1::AbstractSpace{name1,dim1}, space2::AbstractSpace{name2,dim2}) where {name1, name2, dim1, dim2}
        name12 = tuple(name1..., name2...)
        dims = tuple(dim1..., dim2...)
        new{name12, dims}()
    end
end

"""
    FockSpace(name::Symbol)

Create an infinite-dimensional Fock space (number state basis).
Equivalent to `HilbertSpace(name, nothing)`.

# Examples
```jldoctest
julia> F = FockSpace(:F)
ℋ(F)
```
"""
FockSpace(A::Symbol) = HilbertSpace(A, nothing)

"""
    ⊗(space1::AbstractSpace, space2::AbstractSpace)

Tensor product of two Hilbert spaces.
"""
⊗(space1::AbstractSpace, space2::AbstractSpace) = CompositeSpace(space1, space2)

@eval begin
    function Base.$(:(==))(A::AbstractSpace{nameA, dimA}, B::AbstractSpace{nameB, dimB}) where {nameA,nameB,dimA,dimB}
        nameA == nameB && dimA == dimB
    end
end
