# Tensor products: ⊗ operator for kets

# Note: ⊗ is exported from spaces.jl

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
