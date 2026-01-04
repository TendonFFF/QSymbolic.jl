# Inner products: Bra * Ket contractions

# Exports
export InnerProduct

# ==================== INNER PRODUCTS ====================
# Bra * Ket contraction rules

# InnerProduct type for symbolic representation
struct InnerProduct{B1<:AbstractBasis, B2<:AbstractBasis}
    bra::Bra{B1}
    ket::Ket{B2}
end

function Base.show(io::IO, ip::InnerProduct)
    print(io, "⟨", isnothing(ip.bra.index) ? "ψ" : ip.bra.index, "|", 
          isnothing(ip.ket.index) ? "ψ" : ip.ket.index, "⟩")
end

@eval begin
    # Same basis: use contraction rule (default is orthonormal)
    # Handles arbitrary number of indices
    function Base.$(:(*))(bra::Bra{B}, ket::Ket{B}) where B
        # Use the contraction rule system
        return apply_contraction_rule(B, bra.index, ket.index)
    end
    
    # Cross-basis: return symbolic or compute via transform
    function Base.$(:(*))(bra::Bra{B1}, ket::Ket{B2}) where {B1, B2}
        space(B1) == space(B2) || throw(DimensionMismatch("Bra and ket are in different spaces"))
        # Try to find transform
        if has_transform(B2, B1)
            # Transform ket to bra's basis, then compute
            ket_transformed = transform(ket, B1)
            return bra * ket_transformed
        elseif has_transform(B1, B2)
            # Transform bra to ket's basis (priority to ket)
            bra_transformed = transform(Ket(bra), B2)'
            return bra_transformed * ket
        else
            return InnerProduct(bra, ket)
        end
    end
    
    # Bra * WeightedKet → multiply by weight
    Base.$(:(*))(bra::Bra{B}, wk::WeightedKet{B}) where B = simplify((bra * wk.ket) * wk.weight)
    function Base.$(:(*))(bra::Bra{B1}, wk::WeightedKet{B2}) where {B1, B2}
        result = bra * wk.ket
        return result isa Number || result isa AbstractSymbolic ? simplify(result * wk.weight) : result
    end
    
    # Weighted Bra * Ket → multiply by weight
    Base.$(:(*))(wb::WeightedBra{B}, ket::Ket{B}) where B = simplify((wb.bra * ket) * wb.weight)
    function Base.$(:(*))(wb::WeightedBra{B1}, ket::Ket{B2}) where {B1, B2}
        result = wb.bra * ket
        return result isa Number || result isa AbstractSymbolic ? simplify(result * wb.weight) : result
    end
    
    # WeightedBra * WeightedKet → multiply both weights
    Base.$(:(*))(wb::WeightedBra{B}, wk::WeightedKet{B}) where B = simplify((wb.bra * wk.ket) * wb.weight * wk.weight)
    function Base.$(:(*))(wb::WeightedBra{B1}, wk::WeightedKet{B2}) where {B1, B2}
        result = wb.bra * wk.ket
        return result isa Number || result isa AbstractSymbolic ? simplify(result * wb.weight * wk.weight) : result
    end
    
    # SumBra * SumKet → sum over all pairs
    function Base.$(:(*))(sb::SumBra{B1}, sk::SumKet{B2}) where {B1, B2}
        B1 == B2 || space(B1) == space(B2) || throw(DimensionMismatch("Bra and ket are in different spaces"))
        
        result = 0
        for (bra, bra_w) in zip(sb.bras, sb.weights)
            for (ket, ket_w) in zip(sk.kets, sk.weights)
                inner = bra * ket
                if inner isa InnerProduct
                    # Cannot evaluate symbolically, keep as is
                    return inner
                else
                    result += simplify(inner * bra_w * ket_w)
                end
            end
        end
        return simplify(result)
    end
    
    # Mixed cases: promote to Sum types
    Base.$(:(*))(bra::Bra, sk::SumKet) = SumBra([bra], [1]) * sk
    Base.$(:(*))(sb::SumBra, ket::Ket) = sb * SumKet([ket], [1])
    Base.$(:(*))(bra::ProductBra, sk::SumKet) = SumBra([bra], [1]) * sk
    Base.$(:(*))(sb::SumBra, ket::ProductKet) = sb * SumKet([ket], [1])
    
    Base.$(:(*))(bra::Bra, wk::WeightedKet) = bra * wk.ket * wk.weight
    Base.$(:(*))(wb::WeightedBra, ket::Ket) = wb.bra * ket * wb.weight
    Base.$(:(*))(wb::WeightedBra, sk::SumKet) = SumBra([wb.bra], [wb.weight]) * sk
    Base.$(:(*))(sb::SumBra, wk::WeightedKet) = sb * SumKet([wk.ket], [wk.weight])
    
    # ProductBra * ProductKet (same basis) - factorized inner product
    # Since kets and bras are canonically ordered, element-wise multiplication is correct
    Base.$(:(*))(pb::ProductBra{Bs}, pk::ProductKet{Bs}) where {Bs} = 
        simplify(prod(pb.bras[i] * pk.kets[i] for i in 1:length(pb.bras)))
    
    # ProductBra * ProductKet (cross-basis) - try factorized transform
    # Order-independence: both are canonically ordered, so element-wise matching works
    function Base.$(:(*))(pb::ProductBra{Bs1}, pk::ProductKet{Bs2}) where {Bs1,Bs2}
        space(CompositeBasis{Bs1}) == space(CompositeBasis{Bs2}) || 
            throw(DimensionMismatch("Bra and ket are in different spaces"))
        
        # Try factorized: ⟨a₁|⊗⟨a₂|⊗... × |b₁⟩⊗|b₂⟩⊗... = ⟨a₁|b₁⟩ × ⟨a₂|b₂⟩ × ...
        if length(pb.bras) != length(pk.kets)
            throw(DimensionMismatch("Bra and ket have different number of subsystems"))
        end
        
        result = 1
        for i in 1:length(pb.bras)
            result_i = pb.bras[i] * pk.kets[i]
            result = simplify(result * result_i)
        end
        return result
    end
    
    # ==================== PARTIAL CONTRACTION ====================
    # Bra on single space acting on ProductKet - partial trace/contraction
    # Only applies when bra's space matches ONE subsystem of ProductKet
    
    # Helper: Find matching ket in ProductKet for a given bra's space
    function _find_matching_ket(pk::ProductKet, bra_space)
        for (i, k) in enumerate(pk.kets)
            if space(basis(typeof(k))) == bra_space
                return i
            end
        end
        return nothing
    end
    
    # Bra{B} * ProductKet - dispatch based on space relationship
    function Base.$(:(*))(bra::Bra{B}, pk::ProductKet{Bs}) where {B<:Basis, Bs}
        bra_space = space(B)
        pk_space = space(CompositeBasis{Bs})
        
        # Case 1: Same composite space (hybridized ↔ localized) - use transform
        if bra_space == pk_space
            # Try to find transform from CompositeBasis to B (hybridized)
            CB = CompositeBasis{Bs}
            if has_transform(CB, B)
                # Transform ProductKet to hybridized basis
                pk_transformed = transform(Ket(pk), B)
                return bra * pk_transformed
            elseif has_transform(B, CB)
                # Transform bra to localized basis
                bra_transformed = transform(Ket(bra), CB)'
                return bra_transformed * pk
            else
                throw(ArgumentError("No transform defined between hybridized basis $(B) and localized basis $(CB). Define a transform with define_transform!"))
            end
        end
        
        # Case 2: Bra space matches a subsystem - partial contraction
        idx = _find_matching_ket(pk, bra_space)
        
        if isnothing(idx)
            throw(ArgumentError("Bra space $(bra_space) not found in ProductKet subsystems"))
        end
        
        # Extract the matching ket and the rest
        target_ket = pk.kets[idx]
        other_kets = [pk.kets[i] for i in 1:length(pk.kets) if i != idx]
        
        # Compute inner product with matching ket
        inner = bra * target_ket
        
        if inner isa Number && iszero(inner)
            return 0
        end
        
        # Result: inner * (remaining product state)
        if isempty(other_kets)
            return inner
        elseif length(other_kets) == 1
            return inner * other_kets[1]
        else
            return inner * ProductKet(other_kets)
        end
    end
    
    # ==================== LOCALIZED ↔ HYBRIDIZED BASIS (via transform) ====================
    # ProductBra{localized bases} * Ket{hybridized basis on CompositeSpace}
    # Same space, different basis structure - must use transform
    function Base.$(:(*))(pb::ProductBra{Bs}, ket::Ket{B}) where {Bs, B<:Basis{<:CompositeSpace}}
        CB = CompositeBasis{Bs}
        # Verify spaces match
        space(CB) == space(B) || 
            throw(DimensionMismatch("ProductBra and Ket are in different spaces"))
        
        # Try to find transform
        if has_transform(B, CB)
            # Transform ket to localized basis
            ket_transformed = transform(ket, CB)
            return pb * ket_transformed
        elseif has_transform(CB, B)
            # Transform ProductBra to hybridized basis
            pb_transformed = transform(Ket(pb), B)'
            return pb_transformed * ket
        else
            throw(ArgumentError("No transform defined between hybridized basis $(B) and localized basis $(CB). Define a transform with define_transform!"))
        end
    end
end
