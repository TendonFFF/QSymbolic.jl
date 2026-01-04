module QSymbolic

# Re-export useful Symbolics.jl functionality
using Reexport
@reexport using Symbolics: @variables, @syms

# Exports are now in individual files for better organization and maintainability
# See each file for what it exports

# Include source files in dependency order
# 1. Abstract types and foundational types
include("abstract_types.jl")
include("symbolic_scalar_symbolics.jl")  # Symbolics.jl-based implementation
include("spaces.jl")
include("basis.jl")
include("contraction_rules.jl")  # Contraction rule system

# 2. Type definitions (all structs in single files for proper include order)
include("types/states.jl")        # ALL ket/bra types: Ket, Bra, Weighted*, Sum*, Product*
include("types/operators.jl")     # ALL operator types: Outer, Operator, Identity, FunctionOperator, OperatorSum

# 3. Transforms (needs types to be defined)
include("basis_transforms.jl")

# 4. Arithmetic operations (reorganized)
include("arithmetic/ket_arithmetic.jl")      # Adjoint, scalar mult, addition
include("arithmetic/inner_products.jl")      # Bra * Ket contractions
include("arithmetic/tensor_products.jl")     # ⊗ for kets
include("arithmetic/operator_ket.jl")        # Operator * Ket
include("arithmetic/operator_operator.jl")   # Operator * Operator
include("arithmetic/bra_operator.jl")        # Bra * Operator

# 5. Miscellaneous utilities
include("miscellaneous.jl")

end
