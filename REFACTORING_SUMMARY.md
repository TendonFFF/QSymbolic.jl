# QSymbolic.jl Ket Structure Refactoring

## Summary of Changes

This refactoring implements a complete restructuring of the quantum state type hierarchy in QSymbolic.jl as requested.

### New Type Hierarchy

The new hierarchy follows a cleaner structure:

```
AbstractKet
├── Ket                    # Basic ket with index
├── ProductKet             # Tensor product of two kets  
├── WeightedKet            # Ket or ProductKet multiplied by scalar
└── SumKet                 # Linear combination of kets with weights

AbstractBra
├── Bra                    # Lazy adjoint of Ket
├── ProductBra             # Lazy adjoint of ProductKet
├── WeightedBra            # Lazy adjoint of WeightedKet (complex conjugate weight)
└── SumBra                 # Lazy adjoint of SumKet (complex conjugate weights)
```

### Key Design Principles

1. **Basic Elements**: `Ket` and `ProductKet` are the fundamental building blocks at the same hierarchy level

2. **Automatic Promotion**: 
   - Scalar multiplication promotes to `WeightedKet`
   - Addition/subtraction promotes to `SumKet`
   - Simplification happens at the low level (no nested simplification)

3. **Lazy Bras**: Bras are adjoint representations of kets with the same field structure, weights are complex conjugated

4. **Space and Basis**: All objects carry space and basis information for proper identification

5. **Flexible Indices**: Supports arbitrary number of indices - not limited to single or pairs
   - Can be: `n`, `(n,)`, `(n, m)`, `(n, m, k)`, etc.
   - Contraction rules adapt to multi-index nature

6. **Custom Contraction Rules**: Users can define basis-specific inner product behavior
   - Default: simple index checking for orthonormality
   - Custom rules via `define_contraction_rule!`
   - Returns numbers, symbolic expressions, or Dirac deltas

### Implementation Details

#### Core Files Modified/Created

1. **src/states.jl**: New type definitions for Ket, Bra, ProductKet/Bra, WeightedKet/Bra, SumKet/Bra
   - Flexible index type: `KetIndex = Union{SingleIndexValue, Tuple{Vararg{SingleIndexValue}}}`
   - All constructors support arbitrary tuple lengths
   - Show methods for nice display

2. **src/arithmetic.jl**: Arithmetic operations implementing the promotion rules
   - Scalar multiplication: `Number * Ket → WeightedKet`
   - Addition: `Ket + Ket → SumKet`  
   - Adjoint: Complex conjugation of weights
   - Simplification at each operation level

3. **src/contraction_rules.jl**: New file for custom contraction rule system
   - `define_contraction_rule!(basis, rule)`: Define custom rule for a basis
   - `has_contraction_rule(basis)`: Check if custom rule exists
   - `apply_contraction_rule(basis, bra_idx, ket_idx)`: Apply rule (custom or default)
   - Default rule: orthonormal with symbolic Kronecker deltas
   - Multi-index: component-wise checking

4. **src/operators.jl**, **src/basis_transforms.jl**, **src/miscellaneous.jl**: 
   - Updated to use new type names (Ket instead of the previous naming)

#### Implementation Details

The new implementation provides a clean, consistent type hierarchy focused on flexibility and extensibility.

### Usage Examples

#### Basic Usage

```julia
using QSymbolic

# Create basic kets
H = HilbertSpace(:H, 2)
ψ = Ket(H, :ψ)
ϕ = Ket(H, :ϕ)

# Automatic promotion
sum_state = ψ + ϕ           # → SumKet with 2 kets, weights [1, 1]
weighted = 2 * ψ             # → WeightedKet(ψ, 2)
normalized = (ψ + ϕ) / √2    # → SumKet with scaled weights

# Inner products (default orthonormal)
ψ' * ψ  # → 1
ψ' * ϕ  # → 0

# Symbolic indices
n = Sym(:n)
ket_n = Ket(H, n)
ket_n' * ket_n  # → 1
ket_n' * ψ      # → KroneckerDelta(n, :ψ)
```

#### Multi-Index Support

```julia
# Arbitrary number of indices
H1 = HilbertSpace(:A, 2)
H2 = HilbertSpace(:B, 2)  
H3 = HilbertSpace(:C, 2)
composite = H1 ⊗ H2 ⊗ H3

# 3-index ket
ket_3idx = Ket(composite, (:↑, :↓, :↑))

# Also works with symbolic
n, m, k = Sym(:n), Sym(:m), Sym(:k)
ket_symbolic = Ket(composite, (n, m, k))

# Contraction checks all components
ket_3idx' * ket_3idx  # → 1
ket_3idx' * Ket(composite, (:↑, :↓, :↓))  # → 0 (third index differs)
```

#### Custom Contraction Rules

```julia
# Define a non-orthonormal basis
H = HilbertSpace(:H, 2)
Cb = Basis(H, :coherent)

# Coherent states have overlap: ⟨α|β⟩ = exp(-|α-β|²/2)
define_contraction_rule!(typeof(Cb)) do i, j
    if i isa AbstractSymbolic || j isa AbstractSymbolic
        # Return symbolic expression
        exp(-abs(i - j)^2 / 2)
    else
        # Evaluate numerically
        exp(-abs(complex(i) - complex(j))^2 / 2)
    end
end

# Now inner products use custom rule
α = Ket(Cb, 0.5)
β = Ket(Cb, 1.0)
overlap = α' * β  # Uses custom rule

# Clear all custom rules to revert to default
clear_contraction_rules!()
```

#### Basis Transforms

```julia
# Define transformation between bases
H = HilbertSpace(:spin, 2)
Zb = Basis(H, :z)
Xb = Basis(H, :x)

up_z = Ket(Zb, :↑)
down_z = Ket(Zb, :↓)

# Define how x-basis states expand in z-basis
define_transform!(Xb, Zb) do idx
    if idx == :↑
        (up_z + down_z) / √2
    else
        (up_z - down_z) / √2
    end
end

# Cross-basis inner products now work
up_x = Ket(Xb, :↑)
up_z' * up_x  # → 1/√2 (uses transform)
```

### Implementation Details

All types use a consistent, clean interface focused on the core concept that only `Ket` is a true atomic element.

### Testing Status

- [x] Core types compile
- [x] Basic constructors work
- [x] Arithmetic operations implemented
- [x] Contraction rule system implemented
- [ ] Full test suite needs updating for new types
- [ ] Operators need additional testing
- [ ] Basis transforms need testing with new types

### Next Steps

1. Update test suite to use new type names
2. Run full tests to identify any remaining issues
3. Update documentation with new examples
4. Performance testing and optimization
5. Add more examples for custom contraction rules

### Benefits of New Structure

1. **Clearer Hierarchy**: Ket and ProductKet at same level makes sense conceptually
2. **Flexible Indices**: No artificial limit on number of indices
3. **Custom Rules**: Users can define physics-appropriate inner products
4. **Simplified Code**: Promotion rules are explicit and consistent
5. **Better Separation**: Low-level operations on basic types, high-level types expand to basics
6. **Type Safety**: Julia's type system helps catch errors at compile time
