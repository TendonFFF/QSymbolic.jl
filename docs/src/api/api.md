# Full API

This page provides a comprehensive list of all public and documented internal functions in QSymbolic.jl.

For topic-specific documentation, see the individual API pages:
- [Spaces](spaces.md) - Hilbert spaces and composite spaces
- [Bases](bases.md) - Basis types and operations
- [States](states.md) - Kets, bras, and state operations
- [Transforms](transforms.md) - Basis transformations
- [Operators](operators.md) - Quantum operators
- [Symbolic Scalars](symbolic.md) - Symbolic computation

---

## Additional Operator Types

These operator types are documented here for completeness:

```@docs
OperatorSum
OperatorProduct
tr
```

## Contraction Rules (Internal)

Functions for defining custom inner product behavior:

```@docs
define_contraction_rule!
has_contraction_rule
QSymbolic.get_contraction_rule
apply_contraction_rule
clear_contraction_rules!
```

## Basis Utilities (Internal)

```@docs
QSymbolic.bases
QSymbolic.get_transform
```

## Symbolic Utilities (Internal)

```@docs
SymNum
SymExpr
AbstractSymbolic
is_nonnegative
```

## Iteration Protocol

`HilbertSpace` supports iteration for convenient destructuring:

```@docs
Base.iterate(::HilbertSpace)
```

## Tensor Product (Basis)

```@docs
QSymbolic.:⊗(::B1, ::B2) where {B1<:AbstractBasis, B2<:AbstractBasis}
```

## Multiplication

Operator multiplication via bra-ket contraction:

```@docs
Base.:*(::Outer{S}, ::Outer{S}) where S
```

## Index

```@index
```
