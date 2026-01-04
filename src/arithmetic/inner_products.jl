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
end
