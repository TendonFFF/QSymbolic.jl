@doc """
Basis transformation registry and functions.
"""

# Exports
export define_transform!, has_transform, transform, clear_transforms!

# Global registry: (from_basis, to_basis) => transform_function
const BASIS_TRANSFORMS = Dict{Tuple{Type,Type}, Function}()

@doc """
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
up_z, down_z = Ket(Zb, :↑), Ket(Zb, :↓)

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

@doc """
    has_transform(B1, B2) -> Bool

Check if a transform from basis B1 to B2 exists (explicit or factorizable).
"""
has_transform(::Type{B1}, ::Type{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis} = 
    haskey(BASIS_TRANSFORMS, (B1, B2))

# For composite bases: check if factorized transform is possible
function has_transform(::Type{CompositeBasis{As}}, ::Type{CompositeBasis{Bs}}) where {As<:Tuple, Bs<:Tuple}
    # Explicit transform takes priority
    haskey(BASIS_TRANSFORMS, (CompositeBasis{As}, CompositeBasis{Bs})) && return true
    # Must have same number of components
    As_params = As.parameters
    Bs_params = Bs.parameters
    length(As_params) == length(Bs_params) || return false
    # Check if all components can transform (or are identical)
    for (A, B) in zip(As_params, Bs_params)
        (A == B) || has_transform(A, B) || return false
    end
    return true
end

@doc """
    get_transform(B1, B2) -> Function

Get the registered transform function from B1 to B2.
"""
get_transform(::Type{B1}, ::Type{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis} = 
    BASIS_TRANSFORMS[(B1, B2)]

@doc """
    transform(ket, target_basis)

Transform a ket from its current basis to the target basis.
Works for:
- `Ket{B1}` → `B2` (single-system bases, returns SumKet)
- `ProductKet{A1,A2}` → `CompositeBasis{B1,B2}` (factorized transforms)
- `ProductKet{A1,A2}` → `Basis` (composite to eigenbasis)
- `SumKet`, `SumKet`, `WeightedKet` (applies to each component)
"""
function transform(ket::Ket{B1}, ::Type{B2}) where {B1<:AbstractBasis, B2<:AbstractBasis}
    has_transform(B1, B2) || throw(ArgumentError("No transform registered from $B1 to $B2"))
    f = get_transform(B1, B2)
    # Transform function receives the index
    result = f(ket.index)
    result isa AbstractKet || throw(ArgumentError("Transform must return a ket"))
    # Wrap in appropriate sum type if needed
    if result isa Ket
        SumKet(result)
    elseif result isa ProductKet
        SumKet([result], [1])
    else
        result
    end
end

# transform(ProductKet, CompositeBasis) - additional method
function transform(ket::ProductKet{As}, ::Type{CompositeBasis{Bs}}) where {As<:Tuple, Bs<:Tuple}
    As_params = As.parameters
    Bs_params = Bs.parameters
    
    # Identity transform - already in target basis
    if As == Bs
        return SumKet([ket], [1])
    end
    
    # Check for explicit composite transform first
    if haskey(BASIS_TRANSFORMS, (CompositeBasis{As}, CompositeBasis{Bs}))
        f = get_transform(CompositeBasis{As}, CompositeBasis{Bs})
        return f(ket)
    end
    
    # Must have same number of components for factorized transform
    length(As_params) == length(Bs_params) || throw(ArgumentError("Cannot transform between composite bases of different sizes"))
    
    # Factorized transform: (U₁ ⊗ U₂ ⊗ ...) |ψ₁⟩|ψ₂⟩... = (U₁|ψ₁⟩) ⊗ (U₂|ψ₂⟩) ⊗ ...
    transformed_kets = Vector{SumKet}()
    for (i, (A, B)) in enumerate(zip(As_params, Bs_params))
        ki = ket.kets[i]
        if A == B
            push!(transformed_kets, SumKet(ki))
        else
            push!(transformed_kets, transform(ki, B))
        end
    end
    
    # Combine all transformed kets via tensor product expansion
    # Start with first transformed ket
    result_kets = [k for k in transformed_kets[1].kets]
    result_weights = copy(transformed_kets[1].weights)
    
    # Iterate through remaining transformed kets and expand
    for tk in transformed_kets[2:end]
        new_kets = []
        new_weights = promote_type(eltype(result_weights), eltype(tk.weights))[]
        for (k1, w1) in zip(result_kets, result_weights)
            for (k2, w2) in zip(tk.kets, tk.weights)
                # Combine kets into ProductKet
                if k1 isa ProductKet
                    push!(new_kets, ProductKet(vcat(k1.kets, [k2])))
                else
                    push!(new_kets, ProductKet([k1, k2]))
                end
                push!(new_weights, w1 * w2)
            end
        end
        result_kets = new_kets
        result_weights = new_weights
    end
    
    return SumKet(result_kets, result_weights)
end

# transform(ProductKet, Basis) - additional method (for composite space eigenbasis)
function transform(ket::ProductKet{As}, ::Type{B}) where {As<:Tuple, B<:Basis}
    CB = CompositeBasis{As}
    has_transform(CB, B) || throw(ArgumentError("No transform registered from $CB to $B"))
    f = get_transform(CB, B)
    # For ProductKet, create a combined index (tuple of indices from each component)
    combined_idx = Tuple(k.index for k in ket.kets)
    result = f(combined_idx)
    result isa AbstractKet || throw(ArgumentError("Transform must return a ket"))
    if result isa Ket
        SumKet(result)
    else
        result
    end
end

# Transform for SumKet (with CompositeBasis)
function transform(sk::SumKet{CompositeBasis{As},T}, ::Type{B}) where {As<:Tuple,T,B<:AbstractBasis}
    # Identity transform - already in target basis
    if B == CompositeBasis{As}
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

# Transform for SumKet (with simple basis)
function transform(sk::SumKet{B1,T}, ::Type{B2}) where {B1<:Basis,T,B2<:AbstractBasis}
    # Identity transform
    B1 == B2 && return sk
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

# Transform for WeightedKet
function transform(wk::WeightedKet{B1}, ::Type{B2}) where {B1,B2<:AbstractBasis}
    wk.weight * transform(wk.ket, B2)
end

@doc """
    clear_transforms!()

Clear all registered basis transforms.
"""
clear_transforms!() = empty!(BASIS_TRANSFORMS)
