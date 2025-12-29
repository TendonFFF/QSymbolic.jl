@doc "Check that two kets/bras are in the same basis." check_basis
check_basis(::AbstractKet{B1}, ::AbstractKet{B2}) where {B1,B2} =
    B1 == B2 || throw(DimensionMismatch("Kets are in different bases"))
check_basis(::AbstractBra{B1}, ::AbstractBra{B2}) where {B1,B2} =
    B1 == B2 || throw(DimensionMismatch("Bras are in different bases"))

@doc "Check that kets/bras are in the same underlying space." check_space
check_space(k1::AbstractKet, k2::AbstractKet) = 
    space(k1) == space(k2) || throw(DimensionMismatch("Kets are in different spaces"))
check_space(b1::AbstractBra, b2::AbstractBra) = 
    space(b1) == space(b2) || throw(DimensionMismatch("Bras are in different spaces"))

@eval begin
    Base.$(:(==))(ket1::BasisKet{B1}, ket2::BasisKet{B2}) where {B1, B2} =
        B1 == B2 && ket1.index == ket2.index
    Base.$(:(==))(bra1::BasisBra{B1}, bra2::BasisBra{B2}) where {B1, B2} =
        B1 == B2 && bra1.index == bra2.index

    function Base.$(:(+))(ket1::weightedKet{B1}, ket2::weightedKet{B2}) where {B1, B2}
        check_basis(ket1, ket2)
        if ket1.Ket == ket2.Ket
            w = ket1.weight + ket2.weight
            return weightedKet(ket1.Ket, w)
        else
            sumKet(BasisKet{B1}[ket1.Ket, ket2.Ket], [ket1.weight, ket2.weight])
        end
    end
    function Base.$(:(+))(bra1::weightedBra{B1}, bra2::weightedBra{B2}) where {B1, B2}
        check_basis(bra1, bra2)
        if bra1.Bra == bra2.Bra
            w = bra1.weight + bra2.weight
            return weightedBra(bra1.Bra, w)
        else
            sumBra(BasisBra{B1}[bra1.Bra, bra2.Bra], [bra1.weight, bra2.weight])
        end
    end

    Base.$(:(+))(Ket1::BasisKet,    ket2::weightedKet) = weightedKet(Ket1, 1) + ket2
    Base.$(:(+))(Ket1::BasisKet,    Ket2::BasisKet)    = weightedKet(Ket1, 1) + weightedKet(Ket2, 1)
    Base.$(:(+))(ket1::weightedKet, Ket2::BasisKet)    = ket1 + weightedKet(Ket2, 1)

    Base.$(:(+))(Bra1::BasisBra,    bra2::weightedBra) = weightedBra(Bra1, 1) + bra2
    Base.$(:(+))(Bra1::BasisBra,    Bra2::BasisBra)    = weightedBra(Bra1, 1) + weightedBra(Bra2, 1)
    Base.$(:(+))(bra1::weightedBra, Bra2::BasisBra)    = bra1 + weightedBra(Bra2, 1)

    Base.$(:(-))(Ket1::BasisKet,    ket2::weightedKet) = weightedKet(Ket1, 1) + (-1 * ket2)
    Base.$(:(-))(Ket1::BasisKet,    Ket2::BasisKet)    = weightedKet(Ket1, 1) + weightedKet(Ket2, -1)
    Base.$(:(-))(ket1::weightedKet, Ket2::BasisKet)    = ket1 + weightedKet(Ket2, -1)

    Base.$(:(-))(Bra1::BasisBra,    bra2::weightedBra) = weightedBra(Bra1, 1) + (-1 * bra2)
    Base.$(:(-))(Bra1::BasisBra,    Bra2::BasisBra)    = weightedBra(Bra1, 1) + weightedBra(Bra2, -1)
    Base.$(:(-))(bra1::weightedBra, Bra2::BasisBra)    = bra1 + weightedBra(Bra2, -1)
end

@eval begin
    Base.$(:(*))(W::Number, ket::BasisKet{B}) where B = weightedKet(ket, W)
    Base.$(:(*))(W::Number, bra::BasisBra{B}) where B = weightedBra(bra, W)
    Base.$(:(*))(ket::BasisKet{B}, W::Number) where B = weightedKet(ket, W)
    Base.$(:(*))(bra::BasisBra{B}, W::Number) where B = weightedBra(bra, W)

    Base.$(:(*))(W::Number, ket::weightedKet{B}) where B = weightedKet(ket.Ket, W * ket.weight)
    Base.$(:(*))(ket::weightedKet{B}, W::Number) where B = weightedKet(ket.Ket, W * ket.weight)
    Base.$(:(*))(W::Number, bra::weightedBra{B}) where B = weightedBra(bra.Bra, W * bra.weight)
    Base.$(:(*))(bra::weightedBra{B}, W::Number) where B = weightedBra(bra.Bra, W * bra.weight)

    Base.$(:(*))(W::Number, ket::sumKet{B,T}) where {B,T} = sumKet(ket.kets, W .* ket.weights; name=ket.display_name)
    Base.$(:(*))(ket::sumKet{B,T}, W::Number) where {B,T} = W * ket
    Base.$(:(*))(W::Number, bra::sumBra{B,T}) where {B,T} = sumBra(bra.bras, W .* bra.weights; name=bra.display_name)
    Base.$(:(*))(bra::sumBra{B,T}, W::Number) where {B,T} = W * bra

    Base.$(:(//))(ket::AbstractKet, W::Number) = ket * (1 // W)
    Base.$(:(//))(bra::AbstractBra, W::Number) = bra * (1 // W)
    Base.$(:(/))(ket::AbstractKet, W::Number) = ket * (1 / W)
    Base.$(:(/))(bra::AbstractBra, W::Number) = bra * (1 / W)
end

@doc """
    FockKet(space::HilbertSpace, n::Int)

Create a Fock state |n⟩ in the given infinite-dimensional Hilbert space.
The space must have `nothing` as its dimension (created via `HilbertSpace(:name)` or `FockSpace(:name)`).

# Examples
```jldoctest
julia> F = FockSpace(:F);

julia> n0 = FockKet(F, 0)  # ground state
|0⟩

julia> n1 = FockKet(F, 1)  # first excited
|1⟩
```
""" FockKet
function FockKet(space::HilbertSpace{T,dim}, n::Int) where {T,dim}
    dim isa Tuple{Nothing} || throw(ArgumentError("Not a valid Fock space (dimension is limited)"))
    BasisKet(space, n)
end

@doc """
    FockBra(space::HilbertSpace, n::Int)

Create a Fock bra ⟨n| in the given infinite-dimensional Hilbert space.

# Examples
```jldoctest
julia> F = FockSpace(:F);

julia> FockBra(F, 0)
⟨0|
```
""" FockBra
function FockBra(space::HilbertSpace{T,dim}, n::Int) where {T,dim}
    dim isa Tuple{Nothing} || throw(ArgumentError("Not a valid Fock space (dimension is limited)"))
    BasisBra(space, n)
end

BasisKet(bra::BasisBra{B}) where B<:AbstractBasis = BasisKet{B}(bra.index)
BasisBra(ket::BasisKet{B}) where B<:AbstractBasis = BasisBra{B}(ket.index)

Base.adjoint(ket::BasisKet) = BasisBra(ket)
Base.adjoint(bra::BasisBra) = BasisKet(bra)

weightedKet(bra::weightedBra{B}) where B = weightedKet(BasisKet(bra.Bra), bra.weight')
weightedBra(ket::weightedKet{B}) where B = weightedBra(BasisBra(ket.Ket), ket.weight')

Base.adjoint(ket::weightedKet) = weightedBra(ket)
Base.adjoint(bra::weightedBra) = weightedKet(bra)

Base.adjoint(ket::sumKet{B}) where B = sumBra(adjoint.(ket.kets), adjoint.(ket.weights); name=ket.display_name)
Base.adjoint(bra::sumBra{B}) where B = sumKet(adjoint.(bra.bras), adjoint.(bra.weights); name=bra.display_name)

# Inner products - same basis (orthonormal)
@eval begin
    # Same basis: orthonormal → δᵢⱼ
    function Base.$(:(*))(bra::BasisBra{B}, ket::BasisKet{B}) where B
        i, j = bra.index, ket.index
        # If either index is symbolic, return Kronecker delta
        if i isa AbstractSymbolic || j isa AbstractSymbolic
            # Try to simplify if indices are equal symbols
            if i == j
                return 1
            else
                return KroneckerDelta(i, j)
            end
        else
            # Both concrete: evaluate directly
            return i == j ? 1 : 0
        end
    end
    
    # Cross-basis: return symbolic or compute via transform
    function Base.$(:(*))(bra::BasisBra{B1}, ket::BasisKet{B2}) where {B1, B2}
        space(B1) == space(B2) || throw(DimensionMismatch("Bra and ket are in different spaces"))
        # Try to find transform
        if has_transform(B2, B1)
            # Transform ket to bra's basis, then compute
            ket_transformed = transform(ket, B1)
            return bra * ket_transformed
        elseif has_transform(B1, B2)
            # Transform bra to ket's basis
            bra_transformed = transform(BasisKet(bra), B2)'
            return bra_transformed * ket
        else
            return InnerProduct(bra, ket)
        end
    end
end

@eval begin
    function Base.$(:(*))(bras::sumBra{B1}, kets::sumKet{B2}) where {B1, B2}
        B1 == B2 || space(B1) == space(B2) || throw(DimensionMismatch("Bra and ket are in different spaces"))
        bra_list = bras.bras
        bra_weights = bras.weights
        ket_list = kets.kets
        ket_weights = kets.weights

        iter = Iterators.product(1:length(bra_list), 1:length(ket_list))

        sum(iter) do (i, j)
            bra = bra_list[i]
            ket = ket_list[j]
            w = bra_weights[i] * ket_weights[j]
            result = bra * ket
            result isa InnerProduct ? result : result * w
        end
    end
    
    # Specific fallbacks to avoid ambiguity
    Base.$(:(*))(bras::BasisBra, kets::weightedKet) = sumBra(bras) * sumKet(kets)
    Base.$(:(*))(bras::BasisBra, kets::sumKet) = sumBra(bras) * kets
    Base.$(:(*))(bras::weightedBra, kets::BasisKet) = sumBra(bras) * sumKet(kets)
    Base.$(:(*))(bras::weightedBra, kets::weightedKet) = sumBra(bras) * sumKet(kets)
    Base.$(:(*))(bras::weightedBra, kets::sumKet) = sumBra(bras) * kets
    Base.$(:(*))(bras::sumBra, kets::BasisKet) = bras * sumKet(kets)
    Base.$(:(*))(bras::sumBra, kets::weightedKet) = bras * sumKet(kets)
end

# Helper to format index (single or multi)
function _format_index(idx)
    if isnothing(idx)
        "ψ"
    elseif idx isa Tuple
        join(string.(idx), ",")
    else
        string(idx)
    end
end

# Show methods
function Base.show(io::IO, ket::BasisKet)
    name = _format_index(ket.index)
    print(io, "|", name, "⟩")
end

function Base.show(io::IO, bra::BasisBra)
    name = _format_index(bra.index)
    print(io, "⟨", name, "|")
end

function Base.show(io::IO, ket::weightedKet)
    print(io, ket.weight, "·", ket.Ket)
end

function Base.show(io::IO, bra::weightedBra)
    print(io, bra.weight, "·", bra.Bra)
end

function Base.show(io::IO, ket::sumKet)
    if !isnothing(ket.display_name)
        print(io, "|", ket.display_name, "⟩")
    else
        for (i, (k, w)) in enumerate(zip(ket.kets, ket.weights))
            i > 1 && print(io, w >= 0 ? " + " : " - ")
            i == 1 && w < 0 && print(io, "-")
            abs(w) != 1 && print(io, abs(w), "·")
            print(io, k)
        end
    end
end

function Base.show(io::IO, bra::sumBra)
    if !isnothing(bra.display_name)
        print(io, "⟨", bra.display_name, "|")
    else
        for (i, (b, w)) in enumerate(zip(bra.bras, bra.weights))
            i > 1 && print(io, w >= 0 ? " + " : " - ")
            i == 1 && w < 0 && print(io, "-")
            abs(w) != 1 && print(io, abs(w), "·")
            print(io, b)
        end
    end
end

function Base.show(io::IO, ip::InnerProduct)
    print(io, "⟨", isnothing(ip.bra.index) ? "ψ" : ip.bra.index, "|", 
          isnothing(ip.ket.index) ? "ψ" : ip.ket.index, "⟩")
end

function Base.show(io::IO, b::Basis{S,name}) where {S,name}
    print(io, "Basis{", name, "}")
end

function Base.show(io::IO, ::CompositeBasis{B1,B2}) where {B1,B2}
    show(io, B1())
    print(io, "⊗")
    show(io, B2())
end

function Base.show(io::IO, pk::ProductKet)
    print(io, pk.ket1, "⊗", pk.ket2)
end

function Base.show(io::IO, pb::ProductBra)
    print(io, pb.bra1, "⊗", pb.bra2)
end

function Base.show(io::IO, sk::SumProductKet)
    if !isnothing(sk.display_name)
        print(io, "|", sk.display_name, "⟩")
    else
        for (i, (k, w)) in enumerate(zip(sk.kets, sk.weights))
            if i > 1
                # For symbolic weights, always use +
                if w isa Number && !(w isa AbstractSymbolic) && real(w) < 0
                    print(io, " - ")
                    w = -w
                else
                    print(io, " + ")
                end
            elseif w isa Number && !(w isa AbstractSymbolic) && real(w) < 0
                print(io, "-")
                w = -w
            end
            if w != 1
                print(io, "(", w, ")·")
            end
            print(io, k)
        end
    end
end

function Base.show(io::IO, sb::SumProductBra)
    if !isnothing(sb.display_name)
        print(io, "⟨", sb.display_name, "|")
    else
        for (i, (b, w)) in enumerate(zip(sb.bras, sb.weights))
            if i > 1
                # For symbolic weights, always use +
                if w isa Number && !(w isa AbstractSymbolic) && real(w) < 0
                    print(io, " - ")
                    w = -w
                else
                    print(io, " + ")
                end
            elseif w isa Number && !(w isa AbstractSymbolic) && real(w) < 0
                print(io, "-")
                w = -w
            end
            if w != 1
                print(io, "(", w, ")·")
            end
            print(io, b)
        end
    end
end

# Tensor product of kets and bras
⊗(k1::BasisKet{B1}, k2::BasisKet{B2}) where {B1,B2} = ProductKet(k1, k2)
⊗(b1::BasisBra{B1}, b2::BasisBra{B2}) where {B1,B2} = ProductBra(b1, b2)

# Adjoint for product states
Base.adjoint(pk::ProductKet) = ProductBra(adjoint(pk.ket1), adjoint(pk.ket2))
Base.adjoint(pb::ProductBra) = ProductKet(adjoint(pb.bra1), adjoint(pb.bra2))
Base.adjoint(sk::SumProductKet) = SumProductBra([adjoint(k) for k in sk.kets], adjoint.(sk.weights); name=sk.display_name)
Base.adjoint(sb::SumProductBra) = SumProductKet([adjoint(b) for b in sb.bras], adjoint.(sb.weights); name=sb.display_name)

# Inner product of product states (same basis)
@eval begin
    Base.$(:(*))(pb::ProductBra{B1,B2}, pk::ProductKet{B1,B2}) where {B1,B2} = 
        (pb.bra1 * pk.ket1) * (pb.bra2 * pk.ket2)
    
    # ProductBra * SumProductKet
    function Base.$(:(*))(pb::ProductBra{B1,B2}, sk::SumProductKet{B1,B2}) where {B1,B2}
        sum(zip(sk.kets, sk.weights)) do (k, w)
            (pb * k) * w
        end
    end
    
    # SumProductBra * ProductKet
    function Base.$(:(*))(sb::SumProductBra{B1,B2}, pk::ProductKet{B1,B2}) where {B1,B2}
        sum(zip(sb.bras, sb.weights)) do (b, w)
            (b * pk) * w
        end
    end
    
    # SumProductBra * SumProductKet
    function Base.$(:(*))(sb::SumProductBra{B1,B2}, sk::SumProductKet{B1,B2}) where {B1,B2}
        iter = Iterators.product(1:length(sb.bras), 1:length(sk.kets))
        sum(iter) do (i, j)
            (sb.bras[i] * sk.kets[j]) * sb.weights[i] * sk.weights[j]
        end
    end
end

# Cross-basis inner products for product states
# When bra and ket are in different bases but same space, transform and compute

# Helper to check if product state bases are compatible (same composite space)
function _product_bases_same_space(::Type{CompositeBasis{A1,A2}}, ::Type{CompositeBasis{B1,B2}}) where {A1,A2,B1,B2}
    space(CompositeBasis{A1,A2}) == space(CompositeBasis{B1,B2})
end

# ProductBra × ProductKet (cross-basis)
function Base.:*(pb::ProductBra{A1,A2}, pk::ProductKet{B1,B2}) where {A1,A2,B1,B2}
    # Check same composite space
    _product_bases_same_space(CompositeBasis{A1,A2}, CompositeBasis{B1,B2}) || 
        throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    # Try to transform ket to bra's basis
    CB_bra = CompositeBasis{A1,A2}
    CB_ket = CompositeBasis{B1,B2}
    
    if has_transform(CB_ket, CB_bra)
        pk_transformed = transform(pk, CB_bra)
        return pb * pk_transformed
    elseif has_transform(CB_bra, CB_ket)
        pb_transformed = transform(pb', CB_ket)'
        return pb_transformed * pk
    else
        # Try factorized: if each subsystem can transform
        # ⟨a₁|⊗⟨a₂| × |b₁⟩⊗|b₂⟩ = ⟨a₁|b₁⟩ × ⟨a₂|b₂⟩
        return (pb.bra1 * pk.ket1) * (pb.bra2 * pk.ket2)
    end
end

# ProductBra × SumProductKet (cross-basis)
function Base.:*(pb::ProductBra{A1,A2}, sk::SumProductKet{B1,B2,T}) where {A1,A2,B1,B2,T}
    _product_bases_same_space(CompositeBasis{A1,A2}, CompositeBasis{B1,B2}) || 
        throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    sum(zip(sk.kets, sk.weights)) do (k, w)
        (pb * k) * w
    end
end

# SumProductBra × ProductKet (cross-basis)
function Base.:*(sb::SumProductBra{A1,A2,T}, pk::ProductKet{B1,B2}) where {A1,A2,T,B1,B2}
    _product_bases_same_space(CompositeBasis{A1,A2}, CompositeBasis{B1,B2}) || 
        throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    sum(zip(sb.bras, sb.weights)) do (b, w)
        (b * pk) * w
    end
end

# SumProductBra × SumProductKet (cross-basis)
function Base.:*(sb::SumProductBra{A1,A2,T1}, sk::SumProductKet{B1,B2,T2}) where {A1,A2,T1,B1,B2,T2}
    _product_bases_same_space(CompositeBasis{A1,A2}, CompositeBasis{B1,B2}) || 
        throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    iter = Iterators.product(1:length(sb.bras), 1:length(sk.kets))
    sum(iter) do (i, j)
        (sb.bras[i] * sk.kets[j]) * sb.weights[i] * sk.weights[j]
    end
end

# Cross-basis: single-space bra × composite ket (eigenbasis ↔ product basis)
# BasisBra{Basis{CompositeSpace}} × ProductKet{B1,B2}
function Base.:*(bra::BasisBra{B}, pk::ProductKet{B1,B2}) where {B<:Basis, B1,B2}
    CB = CompositeBasis{B1,B2}
    space(B) == space(CB) || throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    if has_transform(CB, B)
        pk_transformed = transform(pk, B)
        return bra * pk_transformed
    elseif has_transform(B, CB)
        bra_transformed = transform(BasisKet{B}(bra.index), CB)'
        return bra_transformed * pk
    else
        throw(ArgumentError("No transform between $B and $CB. Define one with define_transform!"))
    end
end

# BasisBra{Basis{CompositeSpace}} × SumProductKet{B1,B2}
function Base.:*(bra::BasisBra{B}, sk::SumProductKet{B1,B2,T}) where {B<:Basis, B1,B2,T}
    CB = CompositeBasis{B1,B2}
    space(B) == space(CB) || throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    if has_transform(CB, B)
        sk_transformed = transform(sk, B)
        return bra * sk_transformed
    elseif has_transform(B, CB)
        bra_transformed = transform(BasisKet{B}(bra.index), CB)'
        return bra_transformed * sk
    else
        throw(ArgumentError("No transform between $B and $CB. Define one with define_transform!"))
    end
end

# sumBra{Basis{CompositeSpace}} × ProductKet{B1,B2}
function Base.:*(bra::sumBra{B}, pk::ProductKet{B1,B2}) where {B<:Basis, B1,B2}
    CB = CompositeBasis{B1,B2}
    space(B) == space(CB) || throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    sum(zip(bra.bras, bra.weights)) do (b, w)
        (b * pk) * w
    end
end

# sumBra{Basis{CompositeSpace}} × SumProductKet{B1,B2}
function Base.:*(bra::sumBra{B}, sk::SumProductKet{B1,B2,T}) where {B<:Basis, B1,B2,T}
    CB = CompositeBasis{B1,B2}
    space(B) == space(CB) || throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    iter = Iterators.product(1:length(bra.bras), 1:length(sk.kets))
    sum(iter) do (i, j)
        (bra.bras[i] * sk.kets[j]) * bra.weights[i] * sk.weights[j]
    end
end

# ProductBra{B1,B2} × BasisKet{Basis{CompositeSpace}}
function Base.:*(pb::ProductBra{B1,B2}, ket::BasisKet{B}) where {B1,B2,B<:Basis}
    CB = CompositeBasis{B1,B2}
    space(B) == space(CB) || throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    if has_transform(B, CB)
        ket_transformed = transform(ket, CB)
        return pb * ket_transformed
    elseif has_transform(CB, B)
        pb_transformed = transform(pb', B)'
        return pb_transformed * ket
    else
        throw(ArgumentError("No transform between $B and $CB. Define one with define_transform!"))
    end
end

# ProductBra{B1,B2} × sumKet{Basis{CompositeSpace}}
function Base.:*(pb::ProductBra{B1,B2}, ket::sumKet{B,T}) where {B1,B2,B<:Basis,T}
    CB = CompositeBasis{B1,B2}
    space(B) == space(CB) || throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    sum(zip(ket.kets, ket.weights)) do (k, w)
        (pb * k) * w
    end
end

# SumProductBra{B1,B2} × BasisKet{Basis{CompositeSpace}}
function Base.:*(sb::SumProductBra{B1,B2,T}, ket::BasisKet{B}) where {B1,B2,T,B<:Basis}
    CB = CompositeBasis{B1,B2}
    space(B) == space(CB) || throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    sum(zip(sb.bras, sb.weights)) do (b, w)
        (b * ket) * w
    end
end

# SumProductBra{B1,B2} × sumKet{Basis{CompositeSpace}}
function Base.:*(sb::SumProductBra{B1,B2,T1}, ket::sumKet{B,T2}) where {B1,B2,T1,B<:Basis,T2}
    CB = CompositeBasis{B1,B2}
    space(B) == space(CB) || throw(DimensionMismatch("Bra and ket are in different spaces"))
    
    iter = Iterators.product(1:length(sb.bras), 1:length(ket.kets))
    sum(iter) do (i, j)
        (sb.bras[i] * ket.kets[j]) * sb.weights[i] * ket.weights[j]
    end
end

# Scalar multiplication for product states
@eval begin
    Base.$(:(*))(W::Number, pk::ProductKet) = SumProductKet([pk], [W])
    Base.$(:(*))(pk::ProductKet, W::Number) = W * pk
    Base.$(:(*))(W::Number, sk::SumProductKet{B1,B2,T}) where {B1,B2,T} = 
        SumProductKet(sk.kets, W .* sk.weights; name=sk.display_name)
    Base.$(:(*))(sk::SumProductKet, W::Number) = W * sk
    Base.$(:(//))(pk::ProductKet, W::Number) = pk * (1 // W)
    Base.$(:(/))(pk::ProductKet, W::Number) = pk * (1 / W)
    Base.$(:(//))(sk::SumProductKet, W::Number) = sk * (1 // W)
    Base.$(:(/))(sk::SumProductKet, W::Number) = sk * (1 / W)
    
    # Scalar multiplication for product bras
    Base.$(:(*))(W::Number, pb::ProductBra) = SumProductBra([pb], [W])
    Base.$(:(*))(pb::ProductBra, W::Number) = W * pb
    Base.$(:(*))(W::Number, sb::SumProductBra{B1,B2,T}) where {B1,B2,T} = 
        SumProductBra(sb.bras, W .* sb.weights; name=sb.display_name)
    Base.$(:(*))(sb::SumProductBra, W::Number) = W * sb
    Base.$(:(//))(pb::ProductBra, W::Number) = pb * (1 // W)
    Base.$(:(/))(pb::ProductBra, W::Number) = pb * (1 / W)
    Base.$(:(//))(sb::SumProductBra, W::Number) = sb * (1 // W)
    Base.$(:(/))(sb::SumProductBra, W::Number) = sb * (1 / W)
end

# Addition/subtraction for product states
@eval begin
    # ProductKet + ProductKet -> SumProductKet
    function Base.$(:(+))(pk1::ProductKet{B1,B2}, pk2::ProductKet{B1,B2}) where {B1,B2}
        SumProductKet([pk1, pk2], [1, 1])
    end
    function Base.$(:(-))(pk1::ProductKet{B1,B2}, pk2::ProductKet{B1,B2}) where {B1,B2}
        SumProductKet([pk1, pk2], [1, -1])
    end
    
    # ProductKet + SumProductKet
    function Base.$(:(+))(pk::ProductKet{B1,B2}, sk::SumProductKet{B1,B2,T}) where {B1,B2,T}
        SumProductKet(vcat([pk], sk.kets), vcat([one(T)], sk.weights))
    end
    function Base.$(:(+))(sk::SumProductKet{B1,B2,T}, pk::ProductKet{B1,B2}) where {B1,B2,T}
        SumProductKet(vcat(sk.kets, [pk]), vcat(sk.weights, [one(T)]))
    end
    function Base.$(:(-))(pk::ProductKet{B1,B2}, sk::SumProductKet{B1,B2,T}) where {B1,B2,T}
        SumProductKet(vcat([pk], sk.kets), vcat([one(T)], -sk.weights))
    end
    function Base.$(:(-))(sk::SumProductKet{B1,B2,T}, pk::ProductKet{B1,B2}) where {B1,B2,T}
        SumProductKet(vcat(sk.kets, [pk]), vcat(sk.weights, [-one(T)]))
    end
    
    # SumProductKet + SumProductKet
    function Base.$(:(+))(sk1::SumProductKet{B1,B2}, sk2::SumProductKet{B1,B2}) where {B1,B2}
        SumProductKet(vcat(sk1.kets, sk2.kets), vcat(sk1.weights, sk2.weights))
    end
    function Base.$(:(-))(sk1::SumProductKet{B1,B2}, sk2::SumProductKet{B1,B2}) where {B1,B2}
        SumProductKet(vcat(sk1.kets, sk2.kets), vcat(sk1.weights, -sk2.weights))
    end
    
    # ProductBra + ProductBra -> SumProductBra
    function Base.$(:(+))(pb1::ProductBra{B1,B2}, pb2::ProductBra{B1,B2}) where {B1,B2}
        SumProductBra([pb1, pb2], [1, 1])
    end
    function Base.$(:(-))(pb1::ProductBra{B1,B2}, pb2::ProductBra{B1,B2}) where {B1,B2}
        SumProductBra([pb1, pb2], [1, -1])
    end
    
    # ProductBra + SumProductBra
    function Base.$(:(+))(pb::ProductBra{B1,B2}, sb::SumProductBra{B1,B2,T}) where {B1,B2,T}
        SumProductBra(vcat([pb], sb.bras), vcat([one(T)], sb.weights))
    end
    function Base.$(:(+))(sb::SumProductBra{B1,B2,T}, pb::ProductBra{B1,B2}) where {B1,B2,T}
        SumProductBra(vcat(sb.bras, [pb]), vcat(sb.weights, [one(T)]))
    end
    function Base.$(:(-))(pb::ProductBra{B1,B2}, sb::SumProductBra{B1,B2,T}) where {B1,B2,T}
        SumProductBra(vcat([pb], sb.bras), vcat([one(T)], -sb.weights))
    end
    function Base.$(:(-))(sb::SumProductBra{B1,B2,T}, pb::ProductBra{B1,B2}) where {B1,B2,T}
        SumProductBra(vcat(sb.bras, [pb]), vcat(sb.weights, [-one(T)]))
    end
    
    # SumProductBra + SumProductBra
    function Base.$(:(+))(sb1::SumProductBra{B1,B2}, sb2::SumProductBra{B1,B2}) where {B1,B2}
        SumProductBra(vcat(sb1.bras, sb2.bras), vcat(sb1.weights, sb2.weights))
    end
    function Base.$(:(-))(sb1::SumProductBra{B1,B2}, sb2::SumProductBra{B1,B2}) where {B1,B2}
        SumProductBra(vcat(sb1.bras, sb2.bras), vcat(sb1.weights, -sb2.weights))
    end
end