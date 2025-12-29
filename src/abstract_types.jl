# Abstract types for the type hierarchy

@doc """
    AbstractSpace{name, dim}

Abstract supertype for all quantum state spaces.

Type parameters:
- `name`: Tuple of symbols identifying the space
- `dim`: Tuple of dimensions (or `nothing` for infinite-dim)
""" AbstractSpace
abstract type AbstractSpace{name, dim} end

@doc """
    AbstractBasis{space}

Abstract supertype for bases defined on a Hilbert space.
""" AbstractBasis
abstract type AbstractBasis{space<:AbstractSpace} end

@doc """
    AbstractKet{basis}

Abstract supertype for ket vectors |ψ⟩ in a specific basis.
""" AbstractKet
abstract type AbstractKet{basis<:AbstractBasis} end

@doc """
    AbstractBra{basis}

Abstract supertype for bra vectors ⟨ψ| in a specific basis.
""" AbstractBra
abstract type AbstractBra{basis<:AbstractBasis} end
