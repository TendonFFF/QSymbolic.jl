module QSymbolic

# Re-export useful Symbolics.jl functionality
using Reexport
@reexport using Symbolics: @variables, @syms

# Types - Spaces
export AbstractSpace, HilbertSpace, CompositeSpace, FockSpace

# Types - Bases
export AbstractBasis, Basis, DefaultBasis, CompositeBasis

# Types - States
export AbstractKet, AbstractBra
export BasisKet, BasisBra, weightedKet, weightedBra, sumKet, sumBra
export ProductKet, ProductBra, SumProductKet, SumProductBra
export InnerProduct
export SingleIndex, MultiIndex, KetIndex  # Index type aliases

# Types - Operators
export AbstractOperator, Operator, SumOperator, ScaledOperator, OperatorProduct
export FunctionOperator, AdjointFunctionOperator
export TensorOperator
export OpKet, OpBra, IdentityOp
export create_ladder_operators
export lift, swap, reorder, partial_trace

# Types - Symbolic scalars (now using Symbolics.jl backend)
export AbstractSymbolic, Sym, SymNum, SymExpr, KroneckerDelta, ScaledDelta
export substitute, evaluate, symbols, is_numeric, simplify
export is_real, is_positive, is_integer, is_nonnegative, assumptions

# Functions - Spaces & Bases
export ⊗, space, basis, basisname, basis1, basis2

# Functions - States
export FockKet, FockBra
export check_space, check_basis

# Functions - Transforms
export define_transform!, has_transform, transform, clear_transforms!

# Include source files in dependency order
include("abstract_types.jl")
include("symbolic_scalar_symbolics.jl")  # New Symbolics.jl-based implementation
include("spaces.jl")
include("basis.jl")
include("states.jl")
include("basis_transforms.jl")
include("arithmetic.jl")
include("operators.jl")
include("miscellaneous.jl")

end