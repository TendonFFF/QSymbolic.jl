# Abstract types for the type hierarchy

# Exports
export AbstractSpace, AbstractBasis, AbstractKet, AbstractBra, AbstractOperator

@doc """
    AbstractSpace{name, dim}

Abstract supertype for all quantum state spaces.

Type parameters:
- `name`: Tuple of symbols identifying the space
- `dim`: Tuple of dimensions (or `nothing` for infinite-dim)
"""
abstract type AbstractSpace{name, dim} end

@doc """
    AbstractBasis{space}

Abstract supertype for bases defined on a Hilbert space.
"""
abstract type AbstractBasis{space<:AbstractSpace} end

@doc """
    AbstractKet{basis}

Abstract supertype for ket vectors |ψ⟩ in a specific basis.
"""
abstract type AbstractKet{basis<:AbstractBasis} end

@doc """
    AbstractBra{basis}

Abstract supertype for bra vectors ⟨ψ| in a specific basis.
"""
abstract type AbstractBra{basis<:AbstractBasis} end

@doc """
    AbstractOperator{S}

Abstract supertype for all quantum operators acting on space `S`.
Operators use bra-ket arithmetic for contraction.

The type parameter `S` is typically an `AbstractSpace` subtype, but can also be 
`Nothing` for space-agnostic lazy containers like `OperatorProduct`.
"""
abstract type AbstractOperator{S} end
