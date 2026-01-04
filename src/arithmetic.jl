# Arithmetic operations for the new ket/bra structure
# All operations promote to WeightedKet/WeightedBra or SumKet/SumBra as needed

# ==================== ADJOINT OPERATIONS ====================

# Ket ↔ Bra conversions
Ket(bra::Bra{B}) where B<:AbstractBasis = Ket{B}(bra.index)
Bra(ket::Ket{B}) where B<:AbstractBasis = Bra{B}(ket.index)

Base.adjoint(ket::Ket) = Bra(ket)
Base.adjoint(bra::Bra) = Ket(bra)

Base.adjoint(pk::ProductKet) = ProductBra([adjoint(k) for k in pk.kets])
Base.adjoint(pb::ProductBra) = ProductKet([adjoint(b) for b in pb.bras])

# WeightedKet ↔ WeightedBra with complex conjugate
WeightedKet(wb::WeightedBra{B}) where B = WeightedKet(Ket(wb.bra), wb.weight')
WeightedBra(wk::WeightedKet{B}) where B = WeightedBra(Bra(wk.ket), wk.weight')

Base.adjoint(wk::WeightedKet) = WeightedBra(wk)
Base.adjoint(wb::WeightedBra) = WeightedKet(wb)

# SumKet ↔ SumBra with complex conjugate weights
Base.adjoint(sk::SumKet{B}) where B = SumBra([adjoint(k) for k in sk.kets], adjoint.(sk.weights); name=sk.display_name)
Base.adjoint(sb::SumBra{B}) where B = SumKet([adjoint(b) for b in sb.bras], adjoint.(sb.weights); name=sb.display_name)

# ==================== SCALAR MULTIPLICATION ====================
# Basic ket * scalar → WeightedKet
# ProductKet * scalar → WeightedKet (wrapping ProductKet)

@eval begin
    # Ket * Number → WeightedKet
    Base.$(:(*))(W::Number, ket::Ket{B}) where B = WeightedKet(ket, W)
    Base.$(:(*))(ket::Ket{B}, W::Number) where B = WeightedKet(ket, W)
    Base.$(:(*))(W::AbstractSymbolic, ket::Ket{B}) where B = WeightedKet(ket, W)
    Base.$(:(*))(ket::Ket{B}, W::AbstractSymbolic) where B = WeightedKet(ket, W)
    
    # ProductKet * Number → WeightedKet
    Base.$(:(*))(W::Number, pk::ProductKet) = WeightedKet(pk, W)
    Base.$(:(*))(pk::ProductKet, W::Number) = WeightedKet(pk, W)
    Base.$(:(*))(W::AbstractSymbolic, pk::ProductKet) = WeightedKet(pk, W)
    Base.$(:(*))(pk::ProductKet, W::AbstractSymbolic) = WeightedKet(pk, W)
    
    # WeightedKet * Number → WeightedKet (multiply weights, simplify)
    Base.$(:(*))(W::Number, wk::WeightedKet{B}) where B = WeightedKet(wk.ket, simplify(W * wk.weight))
    Base.$(:(*))(wk::WeightedKet{B}, W::Number) where B = WeightedKet(wk.ket, simplify(W * wk.weight))
    Base.$(:(*))(W::AbstractSymbolic, wk::WeightedKet{B}) where B = WeightedKet(wk.ket, simplify(W * wk.weight))
    Base.$(:(*))(wk::WeightedKet{B}, W::AbstractSymbolic) where B = WeightedKet(wk.ket, simplify(W * wk.weight))
    
    # SumKet * Number → SumKet (multiply all weights, simplify)
    Base.$(:(*))(W::Number, sk::SumKet{B,T}) where {B,T} = SumKet(sk.kets, [simplify(W * w) for w in sk.weights]; name=sk.display_name)
    Base.$(:(*))(sk::SumKet{B,T}, W::Number) where {B,T} = W * sk
    Base.$(:(*))(W::AbstractSymbolic, sk::SumKet{B,T}) where {B,T} = SumKet(sk.kets, [simplify(W * w) for w in sk.weights]; name=sk.display_name)
    Base.$(:(*))(sk::SumKet{B,T}, W::AbstractSymbolic) where {B,T} = W * sk
    
    # Same for Bras
    Base.$(:(*))(W::Number, bra::Bra{B}) where B = WeightedBra(bra, W)
    Base.$(:(*))(bra::Bra{B}, W::Number) where B = WeightedBra(bra, W)
    Base.$(:(*))(W::AbstractSymbolic, bra::Bra{B}) where B = WeightedBra(bra, W)
    Base.$(:(*))(bra::Bra{B}, W::AbstractSymbolic) where B = WeightedBra(bra, W)
    
    Base.$(:(*))(W::Number, pb::ProductBra) = WeightedBra(pb, W)
    Base.$(:(*))(pb::ProductBra, W::Number) = WeightedBra(pb, W)
    Base.$(:(*))(W::AbstractSymbolic, pb::ProductBra) = WeightedBra(pb, W)
    Base.$(:(*))(pb::ProductBra, W::AbstractSymbolic) = WeightedBra(pb, W)
    
    Base.$(:(*))(W::Number, wb::WeightedBra{B}) where B = WeightedBra(wb.bra, simplify(W * wb.weight))
    Base.$(:(*))(wb::WeightedBra{B}, W::Number) where B = WeightedBra(wb.bra, simplify(W * wb.weight))
    Base.$(:(*))(W::AbstractSymbolic, wb::WeightedBra{B}) where B = WeightedBra(wb.bra, simplify(W * wb.weight))
    Base.$(:(*))(wb::WeightedBra{B}, W::AbstractSymbolic) where B = WeightedBra(wb.bra, simplify(W * wb.weight))
    
    Base.$(:(*))(W::Number, sb::SumBra{B,T}) where {B,T} = SumBra(sb.bras, [simplify(W * w) for w in sb.weights]; name=sb.display_name)
    Base.$(:(*))(sb::SumBra{B,T}, W::Number) where {B,T} = W * sb
    Base.$(:(*))(W::AbstractSymbolic, sb::SumBra{B,T}) where {B,T} = SumBra(sb.bras, [simplify(W * w) for w in sb.weights]; name=sb.display_name)
    Base.$(:(*))(sb::SumBra{B,T}, W::AbstractSymbolic) where {B,T} = W * sb
    
    # Division
    Base.$(:(//))(ket::AbstractKet, W::Number) = ket * (1 // W)
    Base.$(:(//))(bra::AbstractBra, W::Number) = bra * (1 // W)
    Base.$(:(/))(ket::AbstractKet, W::Number) = ket * (1 / W)
    Base.$(:(/))(bra::AbstractBra, W::Number) = bra * (1 / W)
end

# ==================== ADDITION/SUBTRACTION ====================
# Basic ket + Basic ket → SumKet
# Simplification: combine identical kets

@eval begin
    # Ket + Ket → SumKet
    function Base.$(:(+))(k1::Ket{B}, k2::Ket{B}) where B
        check_basis(k1, k2)
        if k1 == k2
            # Same ket: add weights
            return WeightedKet(k1, 2)
        else
            return SumKet([k1, k2], [1, 1])
        end
    end
    
    function Base.$(:(-))(k1::Ket{B}, k2::Ket{B}) where B
        check_basis(k1, k2)
        if k1 == k2
            # Same ket: result is zero... but we can't return 0 for kets
            # Return weighted ket with weight 0
            return WeightedKet(k1, 0)
        else
            return SumKet([k1, k2], [1, -1])
        end
    end
    
    # Ket + WeightedKet → SumKet
    Base.$(:(+))(k::Ket{B}, wk::WeightedKet{B}) where B = WeightedKet(k, 1) + wk
    Base.$(:(+))(wk::WeightedKet{B}, k::Ket{B}) where B = wk + WeightedKet(k, 1)
    Base.$(:(-))(k::Ket{B}, wk::WeightedKet{B}) where B = WeightedKet(k, 1) - wk
    Base.$(:(-))(wk::WeightedKet{B}, k::Ket{B}) where B = wk - WeightedKet(k, 1)
    
    # WeightedKet + WeightedKet → WeightedKet or SumKet
    function Base.$(:(+))(wk1::WeightedKet{B}, wk2::WeightedKet{B}) where B
        check_basis(wk1, wk2)
        if wk1.ket == wk2.ket
            # Same ket: add weights and simplify
            w_sum = simplify(wk1.weight + wk2.weight)
            return WeightedKet(wk1.ket, w_sum)
        else
            return SumKet([wk1.ket, wk2.ket], [wk1.weight, wk2.weight])
        end
    end
    
    function Base.$(:(-))(wk1::WeightedKet{B}, wk2::WeightedKet{B}) where B
        check_basis(wk1, wk2)
        if wk1.ket == wk2.ket
            # Same ket: subtract weights and simplify
            w_diff = simplify(wk1.weight - wk2.weight)
            return WeightedKet(wk1.ket, w_diff)
        else
            return SumKet([wk1.ket, wk2.ket], [wk1.weight, -wk2.weight])
        end
    end
    
    # Ket + SumKet → SumKet
    Base.$(:(+))(k::Ket{B}, sk::SumKet{B}) where B = SumKet(vcat([k], sk.kets), vcat([1], sk.weights))
    Base.$(:(+))(sk::SumKet{B}, k::Ket{B}) where B = SumKet(vcat(sk.kets, [k]), vcat(sk.weights, [1]))
    Base.$(:(-))(k::Ket{B}, sk::SumKet{B}) where B = SumKet(vcat([k], sk.kets), vcat([1], -sk.weights))
    Base.$(:(-))(sk::SumKet{B}, k::Ket{B}) where B = SumKet(vcat(sk.kets, [k]), vcat(sk.weights, [-1]))
    
    # WeightedKet + SumKet → SumKet
    Base.$(:(+))(wk::WeightedKet{B}, sk::SumKet{B}) where B = SumKet(vcat([wk.ket], sk.kets), vcat([wk.weight], sk.weights))
    Base.$(:(+))(sk::SumKet{B}, wk::WeightedKet{B}) where B = SumKet(vcat(sk.kets, [wk.ket]), vcat(sk.weights, [wk.weight]))
    Base.$(:(-))(wk::WeightedKet{B}, sk::SumKet{B}) where B = SumKet(vcat([wk.ket], sk.kets), vcat([wk.weight], -sk.weights))
    Base.$(:(-))(sk::SumKet{B}, wk::WeightedKet{B}) where B = SumKet(vcat(sk.kets, [wk.ket]), vcat(sk.weights, [-wk.weight]))
    
    # SumKet + SumKet → SumKet
    Base.$(:(+))(sk1::SumKet{B}, sk2::SumKet{B}) where B = SumKet(vcat(sk1.kets, sk2.kets), vcat(sk1.weights, sk2.weights))
    Base.$(:(-))(sk1::SumKet{B}, sk2::SumKet{B}) where B = SumKet(vcat(sk1.kets, sk2.kets), vcat(sk1.weights, -sk2.weights))
    
    # Same operations for Bras
    function Base.$(:(+))(b1::Bra{B}, b2::Bra{B}) where B
        check_basis(b1, b2)
        if b1 == b2
            return WeightedBra(b1, 2)
        else
            return SumBra([b1, b2], [1, 1])
        end
    end
    
    function Base.$(:(-))(b1::Bra{B}, b2::Bra{B}) where B
        check_basis(b1, b2)
        if b1 == b2
            return WeightedBra(b1, 0)
        else
            return SumBra([b1, b2], [1, -1])
        end
    end
    
    Base.$(:(+))(b::Bra{B}, wb::WeightedBra{B}) where B = WeightedBra(b, 1) + wb
    Base.$(:(+))(wb::WeightedBra{B}, b::Bra{B}) where B = wb + WeightedBra(b, 1)
    Base.$(:(-))(b::Bra{B}, wb::WeightedBra{B}) where B = WeightedBra(b, 1) - wb
    Base.$(:(-))(wb::WeightedBra{B}, b::Bra{B}) where B = wb - WeightedBra(b, 1)
    
    function Base.$(:(+))(wb1::WeightedBra{B}, wb2::WeightedBra{B}) where B
        check_basis(wb1, wb2)
        if wb1.bra == wb2.bra
            w_sum = simplify(wb1.weight + wb2.weight)
            return WeightedBra(wb1.bra, w_sum)
        else
            return SumBra([wb1.bra, wb2.bra], [wb1.weight, wb2.weight])
        end
    end
    
    function Base.$(:(-))(wb1::WeightedBra{B}, wb2::WeightedBra{B}) where B
        check_basis(wb1, wb2)
        if wb1.bra == wb2.bra
            w_diff = simplify(wb1.weight - wb2.weight)
            return WeightedBra(wb1.bra, w_diff)
        else
            return SumBra([wb1.bra, wb2.bra], [wb1.weight, -wb2.weight])
        end
    end
    
    Base.$(:(+))(b::Bra{B}, sb::SumBra{B}) where B = SumBra(vcat([b], sb.bras), vcat([1], sb.weights))
    Base.$(:(+))(sb::SumBra{B}, b::Bra{B}) where B = SumBra(vcat(sb.bras, [b]), vcat(sb.weights, [1]))
    Base.$(:(-))(b::Bra{B}, sb::SumBra{B}) where B = SumBra(vcat([b], sb.bras), vcat([1], -sb.weights))
    Base.$(:(-))(sb::SumBra{B}, b::Bra{B}) where B = SumBra(vcat(sb.bras, [b]), vcat(sb.weights, [-1]))
    
    Base.$(:(+))(wb::WeightedBra{B}, sb::SumBra{B}) where B = SumBra(vcat([wb.bra], sb.bras), vcat([wb.weight], sb.weights))
    Base.$(:(+))(sb::SumBra{B}, wb::WeightedBra{B}) where B = SumBra(vcat(sb.bras, [wb.bra]), vcat(sb.weights, [wb.weight]))
    Base.$(:(-))(wb::WeightedBra{B}, sb::SumBra{B}) where B = SumBra(vcat([wb.bra], sb.bras), vcat([wb.weight], -sb.weights))
    Base.$(:(-))(sb::SumBra{B}, wb::WeightedBra{B}) where B = SumBra(vcat(sb.bras, [wb.bra]), vcat(sb.weights, [-wb.weight]))
    
    Base.$(:(+))(sb1::SumBra{B}, sb2::SumBra{B}) where B = SumBra(vcat(sb1.bras, sb2.bras), vcat(sb1.weights, sb2.weights))
    Base.$(:(-))(sb1::SumBra{B}, sb2::SumBra{B}) where B = SumBra(vcat(sb1.bras, sb2.bras), vcat(sb1.weights, -sb2.weights))
    
    # ProductKet operations
    function Base.$(:(+))(pk1::ProductKet{Bs}, pk2::ProductKet{Bs}) where {Bs}
        if pk1 == pk2
            return WeightedKet(pk1, 2)
        else
            return SumKet([pk1, pk2], [1, 1])
        end
    end
    
    function Base.$(:(-))(pk1::ProductKet{Bs}, pk2::ProductKet{Bs}) where {Bs}
        if pk1 == pk2
            return WeightedKet(pk1, 0)
        else
            return SumKet([pk1, pk2], [1, -1])
        end
    end
    
    # ProductBra operations
    function Base.$(:(+))(pb1::ProductBra{Bs}, pb2::ProductBra{Bs}) where {Bs}
        if pb1 == pb2
            return WeightedBra(pb1, 2)
        else
            return SumBra([pb1, pb2], [1, 1])
        end
    end
    
    function Base.$(:(-))(pb1::ProductBra{Bs}, pb2::ProductBra{Bs}) where {Bs}
        if pb1 == pb2
            return WeightedBra(pb1, 0)
        else
            return SumBra([pb1, pb2], [1, -1])
        end
    end
end

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

# ==================== TENSOR PRODUCT ====================
# Ket ⊗ Ket → ProductKet

⊗(k1::Ket{B1}, k2::Ket{B2}) where {B1,B2} = ProductKet([k1, k2])
⊗(b1::Bra{B1}, b2::Bra{B2}) where {B1,B2} = ProductBra([b1, b2])

# ProductKet ⊗ Ket → ProductKet (chaining)
⊗(pk::ProductKet, k::Ket) = ProductKet(vcat(pk.kets, [k]))
⊗(k::Ket, pk::ProductKet) = ProductKet(vcat([k], pk.kets))

# ProductBra ⊗ Bra → ProductBra (chaining)
⊗(pb::ProductBra, b::Bra) = ProductBra(vcat(pb.bras, [b]))
⊗(b::Bra, pb::ProductBra) = ProductBra(vcat([b], pb.bras))

# ProductKet ⊗ ProductKet → ProductKet (concatenation)
⊗(pk1::ProductKet, pk2::ProductKet) = ProductKet(vcat(pk1.kets, pk2.kets))
⊗(pb1::ProductBra, pb2::ProductBra) = ProductBra(vcat(pb1.bras, pb2.bras))
