module QSymbolic

# Types - Spaces
export AbstractSpace, HilbertSpace, CompositeSpace, FockSpace

# Types - Bases
export AbstractBasis, Basis, DefaultBasis, CompositeBasis

# Types - States
export AbstractKet, AbstractBra
export BasisKet, BasisBra, weightedKet, weightedBra, sumKet, sumBra
export ProductKet, ProductBra, SumProductKet, SumProductBra
export InnerProduct

# Types - Operators
export Operator, AdjointOperator, OpKet, OpBra, OperatorProduct

# Types - Symbolic scalars
export AbstractSymbolic, Sym, SymNum, SymExpr
export substitute, evaluate, symbols, is_numeric, simplify

# Functions - Spaces & Bases
export ⊗, space, basis, basisname, basis1, basis2

# Functions - States
export FockKet, FockBra
export check_space, check_basis

# Functions - Transforms
export define_transform!, has_transform, transform, clear_transforms!

# Include source files in dependency order
include("abstract_types.jl")
include("symbolic_scalar.jl")
include("spaces.jl")
include("basis.jl")
include("states.jl")
include("basis_transforms.jl")
include("arithmetic.jl")
include("operators.jl")
include("miscellaneous.jl")

end