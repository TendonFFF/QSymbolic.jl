"""
Basis transformation registry and functions.
"""

# Global registry: (from_basis, to_basis) => transform_function
const BASIS_TRANSFORMS = Dict{Tuple{Type,Type}, Function}()

"""
    define_transform!(from::Basis, to::Basis, f::Function)

Register a basis transformation. `f(index)` returns a sumKet in target basis.

# Example
```julia
H = HilbertSpace(:spin, 2)
Zb, Xb = Basis(H, :z), Basis(H, :x)
up_z, down_z = BasisKet(Zb, :↑), BasisKet(Zb, :↓)

define_transform!(Xb, Zb) do idx
    idx == :↑ ? (up_z + down_z) / √2 : (up_z - down_z) / √2
end
```
"""
function define_transform!(f::Function, from::Basis{S1,n1}, to::Basis{S2,n2}) where {S1,S2,n1,n2}
    S1 == S2 || throw(ArgumentError("Bases must be on the same space"))
    BASIS_TRANSFORMS[(typeof(from), typeof(to))] = f
    nothing
end

# Also allow explicit composite basis transforms (for entangled bases)
function define_transform!(f::Function, from::CompositeBasis{B1,B2}, to::CompositeBasis{B3,B4}) where {B1,B2,B3,B4}
    BASIS_TRANSFORMS[(typeof(from), typeof(to))] = f
    nothing
end

"""
    has_transform(B1, B2) -> Bool

Check if a transform from basis B1 to B2 exists (explicit or factorizable).
"""
has_transform(::Type{B1}, ::Type{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis} = 
    haskey(BASIS_TRANSFORMS, (B1, B2))

# For composite bases: check if factorized transform is possible
function has_transform(::Type{CompositeBasis{A1,A2}}, ::Type{CompositeBasis{B1,B2}}) where {A1,A2,B1,B2}
    # Explicit transform takes priority
    haskey(BASIS_TRANSFORMS, (CompositeBasis{A1,A2}, CompositeBasis{B1,B2})) && return true
    # Otherwise check if both components can transform (or are identical)
    can_transform_1 = (A1 == B1) || has_transform(A1, B1)
    can_transform_2 = (A2 == B2) || has_transform(A2, B2)
    return can_transform_1 && can_transform_2
end

"""
    get_transform(B1, B2) -> Function

Get the registered transform function from B1 to B2.
"""
get_transform(::Type{B1}, ::Type{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis} = 
    BASIS_TRANSFORMS[(B1, B2)]

"""
    transform(ket::BasisKet{B1}, ::Type{B2}) -> sumKet{B2}

Transform a ket from basis B1 to basis B2.
"""
function transform(ket::BasisKet{B1}, ::Type{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    has_transform(B1, B2) || throw(ArgumentError("No transform registered from $B1 to $B2"))
    f = get_transform(B1, B2)
    result = f(ket.index)
    result isa AbstractKet || throw(ArgumentError("Transform must return a ket"))
    result isa sumKet ? result : sumKet(result)
end

"""
    transform(ket::ProductKet, ::Type{CompositeBasis{B1,B2}})

Transform a product ket to a different composite basis using factorized transforms.
"""
function transform(ket::ProductKet{A1,A2}, ::Type{CompositeBasis{B1,B2}}) where {A1,A2,B1,B2}
    # Check for explicit composite transform first
    if haskey(BASIS_TRANSFORMS, (CompositeBasis{A1,A2}, CompositeBasis{B1,B2}))
        f = get_transform(CompositeBasis{A1,A2}, CompositeBasis{B1,B2})
        return f(ket)
    end
    
    # Factorized transform: (U₁ ⊗ U₂)|ψ₁⟩|ψ₂⟩ = (U₁|ψ₁⟩) ⊗ (U₂|ψ₂⟩)
    ket1_transformed = A1 == B1 ? sumKet(ket.ket1) : transform(ket.ket1, B1)
    ket2_transformed = A2 == B2 ? sumKet(ket.ket2) : transform(ket.ket2, B2)
    
    # Combine: Σᵢⱼ wᵢwⱼ |i⟩⊗|j⟩
    result_kets = ProductKet{B1,B2}[]
    result_weights = promote_type(eltype(ket1_transformed.weights), eltype(ket2_transformed.weights))[]
    
    for (k1, w1) in zip(ket1_transformed.kets, ket1_transformed.weights)
        for (k2, w2) in zip(ket2_transformed.kets, ket2_transformed.weights)
            push!(result_kets, ProductKet(k1, k2))
            push!(result_weights, w1 * w2)
        end
    end
    
    return SumProductKet(result_kets, result_weights)
end

"""
    clear_transforms!()

Clear all registered basis transforms.
"""
clear_transforms!() = empty!(BASIS_TRANSFORMS)
