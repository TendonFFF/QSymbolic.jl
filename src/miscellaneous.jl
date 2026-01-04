# Pretty printing for quantum types

# Exports
export TensorOperator, OpKet, OpBra, IdentityOp
export create_ladder_operators
export lift, swap, reorder, partial_trace

# Helper to get space name for display
_space_name(::Type{<:AbstractSpace{name,dim}}) where {name,dim} = join(string.(name), "⊗")

# Spaces
function Base.show(io::IO, ::HilbertSpace{name, dim}) where {name, dim}
    n = name[1]
    d = dim[1]
    if isnothing(d)
        print(io, "ℋ(", n, ")")
    else
        print(io, "ℋ(", n, ", dim=", d, ")")
    end
end

function Base.show(io::IO, ::CompositeSpace{name, dim}) where {name, dim}
    parts = ["ℋ($(name[i]))" for i in eachindex(name)]
    print(io, join(parts, " ⊗ "))
end
