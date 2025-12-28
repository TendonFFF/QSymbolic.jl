# QSymbolic.jl

A Julia package for symbolic quantum mechanics computations.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/QSymbolic.jl")
```

Or for local development:
```julia
] dev /path/to/QSymbolic.jl
```

## Features

- **Hilbert Spaces**: Define finite and infinite-dimensional Hilbert spaces
- **Kets and Bras**: Symbolic representation of quantum states in Dirac notation
- **Arithmetic Operations**: Addition, subtraction, and scalar multiplication of states
- **Inner Products**: Compute inner products between bras and kets
- **Tensor Products**: Compose Hilbert spaces using `⊗`

## Quick Start

```julia
using QSymbolic

# Define a 2-dimensional Hilbert space (qubit)
H = HilbertSpace(:H, 2)

# Create basis kets
ψ = BasisKet(H, :ψ)
ϕ = BasisKet(H, :ϕ)

# Linear combinations
superposition = ψ + ϕ
weighted = (1/√2) * ψ + (1/√2) * ϕ

# Inner products
ψ' * ψ  # Returns 1
ψ' * ϕ  # Returns 0 (orthogonal basis states)

# Fock space for quantum harmonic oscillator
F = FockSpace(:F)
n0 = FockKet(F, 0)  # Ground state |0⟩
n1 = FockKet(F, 1)  # First excited state |1⟩

# Tensor products of spaces
H1 = HilbertSpace(:A, 2)
H2 = HilbertSpace(:B, 2)
H_composite = H1 ⊗ H2
```

## Types

| Type | Description |
|------|-------------|
| `HilbertSpace` | A named Hilbert space with optional dimension |
| `CompositeSpace` | Tensor product of multiple Hilbert spaces |
| `FockSpace` | Infinite-dimensional Fock space |
| `BasisKet` | A basis ket vector `\|ψ⟩` |
| `BasisBra` | A basis bra vector `⟨ψ\|` |
| `weightedKet` | A ket multiplied by a scalar |
| `weightedBra` | A bra multiplied by a scalar |
| `sumKet` | A linear combination of kets |
| `sumBra` | A linear combination of bras |

## License

MIT License - see [LICENSE](LICENSE) for details.
