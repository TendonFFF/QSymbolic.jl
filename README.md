# QSymbolic.jl

[![CI](https://github.com/TendonFFF/QSymbolic.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/TendonFFF/QSymbolic.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://TendonFFF.github.io/QSymbolic.jl)

A Julia package for symbolic quantum mechanics computations with a clean, extensible type hierarchy.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/TendonFFF/QSymbolic.jl")
```

Or for local development:
```julia
] dev /path/to/QSymbolic.jl
```

## Features

- **Hilbert Spaces**: Define finite and infinite-dimensional Hilbert spaces
- **Kets and Bras**: Symbolic representation of quantum states in Dirac notation
- **Arithmetic Operations**: Addition, subtraction, and scalar multiplication with automatic type promotion
- **Inner Products**: Compute inner products with support for cross-basis transforms
- **Tensor Products**: Compose Hilbert spaces and states using `⊗` - supports arbitrary number of subsystems
- **Custom Contraction Rules**: Define basis-specific inner product behavior
- **Multi-index Support**: Kets support arbitrary number of indices
- **Symbolic Computation**: Full integration with Symbolics.jl for symbolic expressions

## Quick Start

```julia
using QSymbolic

# Define a 2-dimensional Hilbert space (qubit)
H = HilbertSpace(:H, 2)
Hb = Basis(H, :default)

# Create basis kets
ψ = Ket(Hb, :ψ)
ϕ = Ket(Hb, :ϕ)

# Linear combinations (automatic promotion to SumKet)
superposition = ψ + ϕ
weighted = (1/√2) * ψ + (1/√2) * ϕ

# Inner products
ψ' * ψ  # Returns 1
ψ' * ϕ  # Returns 0 (orthogonal basis states)

# Fock space for quantum harmonic oscillator
F = FockSpace(:F)
n0 = FockKet(F, 0)  # Ground state |0⟩
n1 = FockKet(F, 1)  # First excited state |1⟩

# Tensor products of spaces and states
H1 = HilbertSpace(:A, 2)
H1b = Basis(H1, :default)
H2 = HilbertSpace(:B, 2)
H2b = Basis(H2, :default)
H3 = HilbertSpace(:C, 2)
H3b = Basis(H3, :default)
H_composite = H1 ⊗ H2 ⊗ H3

# Tensor product of kets (supports arbitrary number)
ψ1 = Ket(H1b, :ψ)
ψ2 = Ket(H2b, :ϕ)
ψ3 = Ket(H3b, :χ)
product_state = ψ1 ⊗ ψ2 ⊗ ψ3  # |ψ⟩⊗|ϕ⟩⊗|χ⟩
```

## Type Hierarchy

The package uses a clean, flattened hierarchy where only `Ket` is a basic element:

```julia
AbstractKet
├── Ket              # Basic ket with index (ONLY basic element)
├── ProductKet       # Container: vector of kets for tensor products
├── WeightedKet      # Wraps Ket or ProductKet with scalar weight
└── SumKet           # Container: vector of kets with weights
```

| Type | Description |
|------|-------------|
| `HilbertSpace` | A named Hilbert space with optional dimension |
| `CompositeSpace` | Tensor product of multiple Hilbert spaces |
| `FockSpace` | Infinite-dimensional Fock space |
| `Ket` | A basis ket vector `|ψ⟩` |
| `Bra` | A basis bra vector `⟨ψ|` (lazy adjoint of Ket) |
| `ProductKet` | Tensor product of kets: `|ψ⟩⊗|ϕ⟩⊗...` |
| `ProductBra` | Tensor product of bras: `⟨ψ|⊗⟨ϕ|⊗...` |
| `WeightedKet` | A ket multiplied by a scalar |
| `WeightedBra` | A bra multiplied by a scalar |
| `SumKet` | A linear combination of kets |
| `SumBra` | A linear combination of bras |

## Advanced Features

### Custom Contraction Rules

Define non-orthonormal bases with custom inner products:

```julia
# Coherent states with overlap ⟨α|β⟩ = exp(-|α-β|²/2)
Cb = Basis(H, :coherent)
define_contraction_rule!(typeof(Cb)) do i, j
    exp(-abs(i - j)^2 / 2)
end

α = Ket(Cb, 0.5)
β = Ket(Cb, 1.0)
overlap = α' * β  # Uses custom rule
```

### Multi-index Support

Kets support arbitrary number of indices:

```julia
# Single, double, triple, or any number of indices
ket1 = Ket(H, :n)
ket2 = Ket(H, (n, m))
ket3 = Ket(H, (n, m, k))
ket5 = Ket(H, (n, m, k, l, p))  # No limits!
```

## License

MIT License - see [LICENSE](LICENSE) for details.
