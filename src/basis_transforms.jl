"""
Basis transformation registry and functions.
"""

# Global registry: (from_basis, to_basis) => transform_function
const BASIS_TRANSFORMS = Dict{Tuple{Type,Type}, Function}()

"""
    define_transform!(f::Function, from::AbstractBasis, to::AbstractBasis)

Register a basis transformation. `f(index)` or `f(ket)` returns a ket in target basis.

The two bases must be on the same space. This works for:
- `Basis` to `Basis` (same space)
- `CompositeBasis` to `CompositeBasis` (same composite space)
- `Basis` on composite space to `CompositeBasis` (eigenbasis to product basis)
- `CompositeBasis` to `Basis` on composite space (product basis to eigenbasis)

# Example
```julia
H = HilbertSpace(:spin, 2)
Zb, Xb = Basis(H, :z), Basis(H, :x)
up_z, down_z = BasisKet(Zb, :↑), BasisKet(Zb, :↓)

define_transform!(Xb, Zb) do idx
    idx == :↑ ? (up_z + down_z) / √2 : (up_z - down_z) / √2
end

# For composite systems, can transform between eigenbasis and product basis:
composite = H1 ⊗ H2
eigenbasis = Basis(composite, :dressed)
product_basis = B1 ⊗ B2

define_transform!(eigenbasis, product_basis) do idx
    # Return superposition in product basis
    ...
end
```
"""
function define_transform!(f::Function, from::B1, to::B2) where {B1<:AbstractBasis, B2<:AbstractBasis}
    space(B1) == space(B2) || throw(ArgumentError("Bases must be on the same space: $(space(B1)) vs $(space(B2))"))
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
    transform(ket::BasisKet{B1}, ::Type{B2}) -> AbstractKet{B2}

Transform a ket from basis B1 to basis B2.
For single-system bases, returns sumKet. For composite systems, may return SumProductKet.
"""
function transform(ket::BasisKet{B1}, ::Type{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    has_transform(B1, B2) || throw(ArgumentError("No transform registered from $B1 to $B2"))
    f = get_transform(B1, B2)
    # Transform function receives the index
    result = f(ket.index)
    result isa AbstractKet || throw(ArgumentError("Transform must return a ket"))
    # Wrap in appropriate sum type if needed
    if result isa BasisKet
        sumKet(result)
    elseif result isa ProductKet
        SumProductKet([result], [1])
    else
        result
    end
end

"""
    transform(ket::ProductKet, ::Type{CompositeBasis{B1,B2}})

Transform a product ket to a different composite basis using factorized transforms.
"""
function transform(ket::ProductKet{A1,A2}, ::Type{CompositeBasis{B1,B2}}) where {A1,A2,B1,B2}
    # Identity transform - already in target basis
    if A1 == B1 && A2 == B2
        return SumProductKet([ket], [1])
    end
    
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
    transform(ket::ProductKet, ::Type{B}) where {B<:Basis}

Transform a product ket to a single basis on the composite space (e.g., eigenbasis).
Requires an explicit transform from the composite basis to the target basis.
"""
function transform(ket::ProductKet{A1,A2}, ::Type{B}) where {A1,A2,B<:Basis}
    CB = CompositeBasis{A1,A2}
    has_transform(CB, B) || throw(ArgumentError("No transform registered from $CB to $B"))
    f = get_transform(CB, B)
    # For ProductKet, create a combined index (tuple of indices from each component)
    combined_idx = (ket.ket1.index, ket.ket2.index)
    result = f(combined_idx)
    result isa AbstractKet || throw(ArgumentError("Transform must return a ket"))
    if result isa BasisKet
        sumKet(result)
    else
        result
    end
end

# Transform for SumProductKet
function transform(sk::SumProductKet{A1,A2,T}, ::Type{B}) where {A1,A2,T,B<:AbstractBasis}
    # Identity transform - already in target basis
    if B == CompositeBasis{A1,A2}
        return sk
    end
    # Transform each component and combine
    total = nothing
    for (k, w) in zip(sk.kets, sk.weights)
        transformed = transform(k, B)
        term = w * transformed
        total = isnothing(total) ? term : total + term
    end
    isnothing(total) ? 0 : total
end

# Transform for sumKet
function transform(sk::sumKet{B1,T}, ::Type{B2}) where {B1,T,B2<:AbstractBasis}
    has_transform(B1, B2) || throw(ArgumentError("No transform registered from $B1 to $B2"))
    # Transform each component and combine
    total = nothing
    for (k, w) in zip(sk.kets, sk.weights)
        transformed = transform(k, B2)
        term = w * transformed
        total = isnothing(total) ? term : total + term
    end
    isnothing(total) ? 0 : total
end

# Transform for weightedKet
function transform(wk::weightedKet{B1}, ::Type{B2}) where {B1,B2<:AbstractBasis}
    wk.weight * transform(wk.Ket, B2)
end

"""
    clear_transforms!()

Clear all registered basis transforms.
"""
clear_transforms!() = empty!(BASIS_TRANSFORMS)
