# Pretty printing for quantum types

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

# Kets
function Base.show(io::IO, ket::BasisKet{space}) where space
    name = isnothing(ket.display_name) ? "ψ" : string(ket.display_name)
    print(io, "|", name, "⟩")
end

function Base.show(io::IO, ket::weightedKet{space}) where space
    name = isnothing(ket.Ket.display_name) ? "ψ" : string(ket.Ket.display_name)
    w = ket.weight
    if w == 1
        print(io, "|", name, "⟩")
    elseif w == -1
        print(io, "-|", name, "⟩")
    else
        print(io, "(", w, ")|", name, "⟩")
    end
end

function Base.show(io::IO, ket::sumKet{space}) where space
    if !isnothing(ket.display_name)
        print(io, "|", ket.display_name, "⟩")
        return
    end
    
    terms = String[]
    for (k, w) in zip(ket.kets, ket.weights)
        name = isnothing(k.display_name) ? "ψ" : string(k.display_name)
        if w == 1
            push!(terms, "|$name⟩")
        elseif w == -1
            push!(terms, "-|$name⟩")
        else
            push!(terms, "($w)|$name⟩")
        end
    end
    
    result = terms[1]
    for t in terms[2:end]
        if startswith(t, "-")
            result *= " - " * t[2:end]
        else
            result *= " + " * t
        end
    end
    print(io, result)
end

# Bras
function Base.show(io::IO, bra::BasisBra{space}) where space
    name = isnothing(bra.display_name) ? "ψ" : string(bra.display_name)
    print(io, "⟨", name, "|")
end

function Base.show(io::IO, bra::weightedBra{space}) where space
    name = isnothing(bra.Bra.display_name) ? "ψ" : string(bra.Bra.display_name)
    w = bra.weight
    if w == 1
        print(io, "⟨", name, "|")
    elseif w == -1
        print(io, "-⟨", name, "|")
    else
        print(io, "(", w, ")⟨", name, "|")
    end
end

function Base.show(io::IO, bra::sumBra{space}) where space
    if !isnothing(bra.display_name)
        print(io, "⟨", bra.display_name, "|")
        return
    end
    
    terms = String[]
    for (b, w) in zip(bra.bras, bra.weights)
        name = isnothing(b.display_name) ? "ψ" : string(b.display_name)
        if w == 1
            push!(terms, "⟨$name|")
        elseif w == -1
            push!(terms, "-⟨$name|")
        else
            push!(terms, "($w)⟨$name|")
        end
    end
    
    result = terms[1]
    for t in terms[2:end]
        if startswith(t, "-")
            result *= " - " * t[2:end]
        else
            result *= " + " * t
        end
    end
    print(io, result)
end
