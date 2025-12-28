# Abstract types for the type hierarchy

"""
    AbstractSpace{name, dim}

Abstract supertype for all quantum state spaces.

Type parameters:
- `name`: Tuple of symbols identifying the space
- `dim`: Tuple of dimensions (or `nothing` for infinite-dim)
"""
abstract type AbstractSpace{name, dim} end

"""
    AbstractBasis{space}

Abstract supertype for bases defined on a Hilbert space.
"""
abstract type AbstractBasis{space<:AbstractSpace} end

"""
    AbstractKet{basis}

Abstract supertype for ket vectors |ψ⟩ in a specific basis.
"""
abstract type AbstractKet{basis<:AbstractBasis} end

"""
    AbstractBra{basis}

Abstract supertype for bra vectors ⟨ψ| in a specific basis.
"""
abstract type AbstractBra{basis<:AbstractBasis} end
