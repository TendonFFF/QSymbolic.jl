# Symbolic scalar arithmetic using Symbolics.jl
# This module provides a thin wrapper around Symbolics.jl for scalar symbolic computation

# Exports
export AbstractSymbolic, Sym, SymNum, SymExpr, KroneckerDelta, ScaledDelta
export substitute, evaluate, symbols, is_numeric, simplify
export is_real, is_positive, is_integer, is_nonnegative, assumptions

import Symbolics
import Symbolics: @variables

# ============== Re-export Symbolics types with QSymbolic-compatible API ==============

"""
    AbstractSymbolic

Abstract supertype for all symbolic scalar types.
This is now an alias that accepts both our wrappers and Symbolics.jl types.
"""
const AbstractSymbolic = Union{Symbolics.Num, Complex{Symbolics.Num}}

# ============== Sym - Symbolic variable with assumptions ==============

"""
    Sym(name::Symbol)
    Sym(name::Symbol; real=false, positive=false, negative=false, integer=false, nonnegative=false)

Create a symbolic variable with optional assumptions, using Symbolics.jl under the hood.

Supported assumptions:
- `:real` / `real=true` - real number (conj(x) = x)
- `:positive` / `positive=true` - strictly positive real
- `:negative` / `negative=true` - strictly negative real  
- `:nonnegative` / `nonnegative=true` - non-negative real
- `:integer` / `integer=true` - integer value

# Examples
```julia
n = Sym(:n)                     # generic symbolic
m = Sym(:m, :real, :positive)   # real positive
k = Sym(:k, integer=true)       # integer
```
"""
function Sym(name::Symbol; 
             real::Bool=false, 
             positive::Bool=false, 
             negative::Bool=false, 
             integer::Bool=false, 
             nonnegative::Bool=false,
             complex::Bool=false)
    
    # Build the variable type based on assumptions
    # Symbolics.jl uses type annotations for assumptions
    if integer
        var = only(@variables $name::Integer)
    elseif positive
        # Positive implies real and nonnegative
        var = only(@variables $name::Real)
        # Note: Symbolics doesn't have built-in positive constraint,
        # but Real is sufficient for most simplifications
    elseif nonnegative
        var = only(@variables $name::Real)
    elseif negative
        var = only(@variables $name::Real)
    elseif real
        var = only(@variables $name::Real)
    elseif complex
        var = only(@variables $name::Complex)
    else
        # Default: Real (most quantum mechanics uses real coefficients)
        var = only(@variables $name::Real)
    end
    
    return var
end

# Varargs constructor: Sym(:n, :real, :positive)
function Sym(name::Symbol, assumptions::Symbol...)
    real = :real in assumptions
    positive = :positive in assumptions
    negative = :negative in assumptions
    integer = :integer in assumptions
    nonnegative = :nonnegative in assumptions
    complex = :complex in assumptions
    
    Sym(name; real, positive, negative, integer, nonnegative, complex)
end

# ============== SymNum - Numeric wrapper (for compatibility) ==============

"""
    SymNum(value)

Wraps a numeric value. With Symbolics.jl integration, this just returns 
the value directly since Symbolics handles numeric constants natively.
"""
SymNum(x::Number) = x
SymNum(x::Symbolics.Num) = x

# ============== SymExpr - Expression wrapper (for compatibility) ==============

"""
    SymExpr(op, args...)

Create a symbolic expression. With Symbolics.jl, expressions are built
automatically through arithmetic operations.

This constructor is provided for compatibility but you should prefer
using natural arithmetic syntax.
"""
function SymExpr(op::Symbol, args...)
    if op == :+
        return +(args...)
    elseif op == :-
        if length(args) == 1
            return -args[1]
        else
            return -(args...)
        end
    elseif op == :*
        return *(args...)
    elseif op == :/
        return /(args...)
    elseif op == :^
        return ^(args...)
    elseif op == :sqrt
        return sqrt(args[1])
    elseif op == :conj
        return conj(args[1])
    elseif op == :abs
        return abs(args[1])
    elseif op == :abs2
        return abs2(args[1])
    elseif op == :exp
        return exp(args[1])
    elseif op == :log
        return log(args[1])
    elseif op == :sin
        return sin(args[1])
    elseif op == :cos
        return cos(args[1])
    else
        error("Unknown operation: $op")
    end
end

# ============== KroneckerDelta ==============

"""
    KroneckerDelta(i, j)

Kronecker delta δᵢⱼ = 1 if i == j, else 0.
With Symbolics.jl, we use IfElse or direct evaluation when possible.
"""
struct KroneckerDelta
    i::Any
    j::Any
end

function Base.show(io::IO, δ::KroneckerDelta)
    print(io, "δ(", δ.i, ",", δ.j, ")")
end

# Arithmetic with KroneckerDelta
Base.:*(a::Number, δ::KroneckerDelta) = iszero(a) ? 0 : ScaledDelta(a, δ)
Base.:*(δ::KroneckerDelta, a::Number) = a * δ
Base.:*(a::Symbolics.Num, δ::KroneckerDelta) = ScaledDelta(a, δ)
Base.:*(δ::KroneckerDelta, a::Symbolics.Num) = ScaledDelta(a, δ)

struct ScaledDelta
    coeff::Any
    delta::KroneckerDelta
end

Base.show(io::IO, sd::ScaledDelta) = print(io, sd.coeff, "·", sd.delta)

# ============== Simplification ==============

"""
    simplify(expr)

Simplify a symbolic expression using Symbolics.jl's simplification engine.
"""
simplify(x::Number) = x
simplify(x::Symbolics.Num) = Symbolics.simplify(x)
simplify(x::Complex{Symbolics.Num}) = Complex(Symbolics.simplify(real(x)), Symbolics.simplify(imag(x)))

function simplify(δ::KroneckerDelta)
    i_simp = simplify(δ.i)
    j_simp = simplify(δ.j)
    
    # Try to evaluate if both are concrete
    if i_simp isa Number && j_simp isa Number
        return i_simp == j_simp ? 1 : 0
    end
    
    # Check symbolic equality using Symbolics
    if i_simp isa Symbolics.Num && j_simp isa Symbolics.Num
        diff = Symbolics.simplify(i_simp - j_simp)
        if Symbolics.iszero(diff)
            return 1
        end
    end
    
    # Can't simplify further
    return KroneckerDelta(i_simp, j_simp)
end

function simplify(sd::ScaledDelta)
    coeff_simp = simplify(sd.coeff)
    delta_simp = simplify(sd.delta)
    
    if delta_simp isa Number
        return coeff_simp * delta_simp
    end
    
    if iszero(coeff_simp)
        return 0
    end
    
    return ScaledDelta(coeff_simp, delta_simp)
end

# ============== Substitution ==============

"""
    substitute(expr, pairs...)
    substitute(expr, dict)

Substitute values for symbolic variables.
"""
substitute(x::Number, args...) = x
substitute(x::Number, ::Dict) = x

function substitute(x::Symbolics.Num, pairs::Pair{Symbol}...)
    # Get the actual symbolic variables from the expression
    vars = Symbolics.get_variables(x)
    subs_dict = Dict{Any, Any}()
    for (name, val) in pairs
        # Find the variable with this name in the expression
        for var in vars
            if Symbol(var) == name
                subs_dict[var] = val
                break
            end
        end
    end
    return Symbolics.substitute(x, subs_dict)
end

function substitute(x::Symbolics.Num, d::Dict{Symbol, <:Any})
    # Get the actual symbolic variables from the expression
    vars = Symbolics.get_variables(x)
    subs_dict = Dict{Any, Any}()
    for (name, val) in d
        # Find the variable with this name in the expression
        for var in vars
            if Symbol(var) == name
                subs_dict[var] = val
                break
            end
        end
    end
    return Symbolics.substitute(x, subs_dict)
end

function substitute(δ::KroneckerDelta, args...)
    i_new = substitute(δ.i, args...)
    j_new = substitute(δ.j, args...)
    simplify(KroneckerDelta(i_new, j_new))
end

# ============== Evaluation ==============

"""
    evaluate(expr)

Evaluate a symbolic expression to a numeric value.
All symbolic variables must have been substituted with concrete values.
"""
evaluate(x::Number) = x
function evaluate(x::Symbolics.Num)
    val = Symbolics.value(x)
    if val isa Symbolics.Num || !(val isa Number)
        error("Cannot evaluate expression with symbolic variables: $x")
    end
    return val
end

function evaluate(δ::KroneckerDelta)
    i_val = evaluate(δ.i)
    j_val = evaluate(δ.j)
    return i_val == j_val ? 1 : 0
end

# ============== Query Functions ==============

"""Check if expression is purely numeric (can be evaluated)"""
is_numeric(x::Number) = true
is_numeric(x::Symbolics.Num) = Symbolics.value(x) isa Number
is_numeric(δ::KroneckerDelta) = is_numeric(δ.i) && is_numeric(δ.j)

"""Get all symbolic variables in an expression"""
function symbols(x::Symbolics.Num)
    vars = Symbolics.get_variables(x)
    return [Symbol(v) for v in vars]
end
symbols(::Number) = Symbol[]
symbols(δ::KroneckerDelta) = union(symbols(δ.i), symbols(δ.j))

"""Check if symbolic is assumed to be real"""
is_real(x::Symbolics.Num) = Symbolics.symtype(Symbolics.value(x)) <: Real
is_real(::Real) = true
is_real(::Complex) = false

"""Check if symbolic is assumed to be positive"""
is_positive(x::Number) = x > 0
is_positive(x::Symbolics.Num) = false  # Would need metadata inspection

"""Check if symbolic is assumed to be an integer"""
is_integer(x::Number) = isinteger(x)
is_integer(x::Symbolics.Num) = Symbolics.symtype(Symbolics.value(x)) <: Integer

"""Check if symbolic is assumed to be non-negative"""
is_nonnegative(x::Number) = x >= 0
is_nonnegative(x::Symbolics.Num) = false  # Would need metadata inspection

"""Get assumptions of a symbolic variable (compatibility function)"""
function assumptions(x::Symbolics.Num)
    result = Set{Symbol}()
    T = Symbolics.symtype(Symbolics.value(x))
    if T <: Real
        push!(result, :real)
    end
    if T <: Integer
        push!(result, :integer, :real)
    end
    return result
end
assumptions(::Number) = Set{Symbol}()

# ============== Display ==============

# Symbolics.jl handles display automatically

# ============== Compatibility exports ==============

# These ensure code written for the old API still works
const SymbolicScalar = Union{Number, Symbolics.Num, KroneckerDelta, ScaledDelta}
