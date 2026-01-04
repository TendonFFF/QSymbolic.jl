# Custom Contraction Rules

QSymbolic.jl allows you to define custom inner product behavior for specific bases. This is essential for modeling non-orthonormal bases, dressed states, or bases with special structure like Jaynes-Cummings states.

## Why Custom Contraction Rules?

By default, QSymbolic.jl assumes orthonormal bases:
- `⟨i|j⟩ = 1` if indices are equal
- `⟨i|j⟩ = 0` if indices differ
- `⟨i|j⟩ = δᵢⱼ` (KroneckerDelta) if either index is symbolic

But many physical situations require different behavior:
- **Coherent states**: `⟨α|β⟩ = exp(-|α-β|²/2)` (non-orthogonal)
- **Dressed states**: Custom overlap depending on quantum numbers
- **Multi-index bases**: Component-wise checking with special rules

## Defining a Custom Rule

Use `define_contraction_rule!` to specify how inner products behave for a basis:

```julia
using QSymbolic

H, Hb = HilbertSpace(:H, 2)
my_basis = Basis(H, :custom)

# Define the contraction rule
define_contraction_rule!(typeof(my_basis)) do bra_idx, ket_idx
    # Your custom logic here
    # Return: number, symbolic expression, or KroneckerDelta
end
```

The rule function receives:
- `bra_idx`: The index from the bra `⟨bra_idx|`
- `ket_idx`: The index from the ket `|ket_idx⟩`

And should return:
- A number (`0`, `1`, `0.5`, etc.)
- A symbolic expression
- A `KroneckerDelta(i, j)` for symbolic orthonormality

## Example: Orthonormal Basis (Default Behavior)

The default rule implements standard orthonormality:

```julia
H, Hb = HilbertSpace(:H, 2)
Zb = Basis(H, :z)

# This is what the default does (no need to define):
define_contraction_rule!(typeof(Zb)) do i, j
    if i isa AbstractSymbolic || j isa AbstractSymbolic
        # Symbolic: return 1 if same symbol, else KroneckerDelta
        return isequal(i, j) ? 1 : KroneckerDelta(i, j)
    else
        # Concrete: direct comparison
        return isequal(i, j) ? 1 : 0
    end
end

up = Ket(Zb, :↑)
down = Ket(Zb, :↓)

up' * up    # → 1
up' * down  # → 0
```

## Example: Non-Orthogonal Coherent States

Coherent states have non-zero overlap:

```julia
H, Hb = HilbertSpace(:oscillator)
Cb = Basis(H, :coherent)

define_contraction_rule!(typeof(Cb)) do α, β
    # ⟨α|β⟩ = exp(-|α-β|²/2) × exp(i×Im(α*conj(β)))
    # Simplified for real coherent amplitudes:
    if α isa AbstractSymbolic || β isa AbstractSymbolic
        # Keep symbolic
        return exp(-abs(α - β)^2 / 2)
    else
        # Evaluate numerically
        return exp(-abs(α - β)^2 / 2)
    end
end

# Coherent states
α = Ket(Cb, 0.5)
β = Ket(Cb, 1.0)

α' * α  # → 1 (self-overlap)
α' * β  # → exp(-0.125) ≈ 0.882 (non-zero overlap!)
```

## Working with Symbolic Indices

When indices are symbolic (using `Sym`), use these patterns:

```julia
# Check if a value is symbolic
x isa AbstractSymbolic  # true for Sym(:x), false for :x or 1

# Compare values safely (works for symbols and numbers)
isequal(a, b)  # Use this instead of ==

# Return KroneckerDelta for symbolic orthonormality
KroneckerDelta(i, j)  # δᵢⱼ - simplifies to 0 or 1 when possible

# Simplify expressions
KroneckerDelta(i, j) |> simplify
```

!!! warning "Use `isequal` not `==`"
    Always use `isequal(a, b)` for comparisons in contraction rules. Using `==` on symbolic expressions can throw errors or return symbolic booleans that break control flow.

## API Reference

### `define_contraction_rule!(basis_type, rule)`

Define a custom contraction rule for a basis type.

```julia
define_contraction_rule!(typeof(my_basis)) do bra_idx, ket_idx
    # Return the inner product result
end
```

### `has_contraction_rule(basis_type)`

Check if a custom rule exists:

```julia
has_contraction_rule(typeof(my_basis))  # → true/false
```

### `clear_contraction_rules!()`

Remove all custom rules, reverting to default orthonormal behavior:

```julia
clear_contraction_rules!()
```

## Best Practices

1. **Always handle both symbolic and concrete cases** - Check `isa AbstractSymbolic` before comparisons

2. **Use `isequal` for comparisons** - Avoids issues with symbolic equality

3. **Simplify results** - Use `|> simplify` to reduce KroneckerDelta expressions

4. **Return consistent types** - Return numbers for concrete cases, symbolic expressions for symbolic cases

5. **Document your rule** - Complex contraction rules benefit from comments explaining the physics

## See Also

- [Getting Started](@ref) - Basic ket/bra operations
- [Composite Systems](@ref) - Tensor products and multi-index states
- [Symbolic Scalars](symbolic.md) - Working with symbolic expressions
