# Symbolic scalar arithmetic - lazy evaluation

@doc """
    AbstractSymbolic

Abstract supertype for all symbolic scalar types.
""" AbstractSymbolic
abstract type AbstractSymbolic <: Number end

@doc """
    Sym(name::Symbol)
    Sym(name::Symbol, assumptions...)
    Sym(name::Symbol; real=false, positive=false, negative=false, integer=false, nonnegative=false, complex=false)

A symbolic variable with optional assumptions. Arithmetic operations build expression trees.

Supported assumptions:
- `:real` - real number (conj(x) = x)
- `:positive` - strictly positive real (implies :real, :nonnegative)
- `:negative` - strictly negative real (implies :real)
- `:nonnegative` - non-negative real (implies :real)
- `:integer` - integer value (implies :real)
- `:complex` - complex number (default if no real assumption)

# Examples
```julia
n = Sym(:n)                     # generic symbolic
m = Sym(:m, :real, :positive)   # real positive
k = Sym(:k, integer=true)       # integer (implies real)

adjoint(n)  # → n*  (conjugate)
adjoint(m)  # → m   (real, so no conjugate)
```
""" Sym
struct Sym <: AbstractSymbolic
    name::Symbol
    assumptions::Set{Symbol}
end

# Helper to expand implied assumptions
function _expand_assumptions(assumptions::Set{Symbol})
    expanded = copy(assumptions)
    if :positive in expanded
        push!(expanded, :real, :nonnegative)
    end
    if :negative in expanded
        push!(expanded, :real)
    end
    if :nonnegative in expanded
        push!(expanded, :real)
    end
    if :integer in expanded
        push!(expanded, :real)
    end
    expanded
end

# Constructors - use vararg with at least one Symbol to avoid method ambiguity
# Basic constructor: Sym(:n) or Sym(:n, :real) or Sym(:n, :real, :positive)
function Sym(name::Symbol, assumptions::Symbol...)
    if isempty(assumptions)
        Sym(name, Set{Symbol}())
    else
        Sym(name, _expand_assumptions(Set{Symbol}(assumptions)))
    end
end

# Keyword constructor: Sym(:n, real=true, positive=true)
function Sym(name::Symbol; 
             real::Bool=false, 
             positive::Bool=false, 
             negative::Bool=false, 
             integer::Bool=false, 
             nonnegative::Bool=false,
             complex::Bool=false)
    assumptions = Set{Symbol}()
    real && push!(assumptions, :real)
    positive && push!(assumptions, :positive)
    negative && push!(assumptions, :negative)
    integer && push!(assumptions, :integer)
    nonnegative && push!(assumptions, :nonnegative)
    complex && push!(assumptions, :complex)
    Sym(name, _expand_assumptions(assumptions))
end

@doc """Get assumptions of a symbolic variable""" assumptions
assumptions(s::Sym) = s.assumptions
assumptions(::AbstractSymbolic) = Set{Symbol}()

# Equality for Sym - two Syms are equal if they have same name and assumptions
Base.:(==)(a::Sym, b::Sym) = a.name == b.name && a.assumptions == b.assumptions
Base.hash(s::Sym, h::UInt) = hash(s.name, hash(s.assumptions, h))

@doc """
    SymNum(value)

Wraps a numeric value in the symbolic system.
""" SymNum
struct SymNum{T<:Number} <: AbstractSymbolic
    value::T
end

# Auto-convert numbers in symbolic context
SymNum(s::AbstractSymbolic) = s

# Assumption query functions (must come after SymNum is defined)
@doc """Check if symbolic is assumed to be real""" is_real
is_real(s::Sym) = :real in s.assumptions
is_real(::SymNum{<:Real}) = true
is_real(::SymNum) = false
is_real(::AbstractSymbolic) = false  # SymExpr - could be improved

@doc """Check if symbolic is assumed to be positive""" is_positive
is_positive(s::Sym) = :positive in s.assumptions
is_positive(s::SymNum{<:Real}) = s.value > 0
is_positive(::AbstractSymbolic) = false

@doc """Check if symbolic is assumed to be an integer""" is_integer
is_integer(s::Sym) = :integer in s.assumptions
is_integer(::SymNum{<:Integer}) = true
is_integer(::AbstractSymbolic) = false

@doc """Check if symbolic is assumed to be non-negative""" is_nonnegative
is_nonnegative(s::Sym) = :nonnegative in s.assumptions
is_nonnegative(s::SymNum{<:Real}) = s.value >= 0
is_nonnegative(::AbstractSymbolic) = false

@doc """
    SymExpr(op, args...)

A symbolic expression representing `op(args...)`.
Operations are stored as symbols: :+, :-, :*, :/, :^, :sqrt, :conj, etc.
""" SymExpr
struct SymExpr <: AbstractSymbolic
    op::Symbol
    args::Tuple{Vararg{AbstractSymbolic}}
    
    function SymExpr(op::Symbol, args...)
        # Wrap any numbers
        wrapped = map(a -> a isa AbstractSymbolic ? a : SymNum(a), args)
        new(op, wrapped)
    end
end

# Display
function Base.show(io::IO, s::Sym)
    print(io, s.name)
    if get(io, :show_assumptions, false) && !isempty(s.assumptions)
        print(io, "{", join(sort(collect(s.assumptions)), ","), "}")
    end
end
Base.show(io::IO, s::SymNum) = print(io, s.value)

function Base.show(io::IO, e::SymExpr)
    op = e.op
    args = e.args
    
    if op == :sqrt
        print(io, "√", "(", args[1], ")")
    elseif op == :conj
        print(io, args[1], "*")
    elseif op == :neg
        print(io, "-", args[1])
    elseif op == :+ 
        print(io, "(", args[1])
        for a in args[2:end]
            print(io, " + ", a)
        end
        print(io, ")")
    elseif op == :-
        print(io, "(", args[1], " - ", args[2], ")")
    elseif op == :*
        for (i, a) in enumerate(args)
            i > 1 && print(io, "·")
            print(io, a)
        end
    elseif op == :/
        print(io, "(", args[1], "/", args[2], ")")
    elseif op == :^
        print(io, args[1], "^", args[2])
    else
        print(io, op, "(", join(args, ", "), ")")
    end
end

# ============== Arithmetic operations ==============

# Addition
Base.:+(a::AbstractSymbolic, b::AbstractSymbolic) = SymExpr(:+, a, b)
Base.:+(a::AbstractSymbolic, b::Number) = SymExpr(:+, a, SymNum(b))
Base.:+(a::Number, b::AbstractSymbolic) = SymExpr(:+, SymNum(a), b)

# Subtraction
Base.:-(a::AbstractSymbolic, b::AbstractSymbolic) = SymExpr(:-, a, b)
Base.:-(a::AbstractSymbolic, b::Number) = SymExpr(:-, a, SymNum(b))
Base.:-(a::Number, b::AbstractSymbolic) = SymExpr(:-, SymNum(a), b)
Base.:-(a::AbstractSymbolic) = SymExpr(:neg, a)

# Multiplication
Base.:*(a::AbstractSymbolic, b::AbstractSymbolic) = SymExpr(:*, a, b)
Base.:*(a::AbstractSymbolic, b::Number) = SymExpr(:*, a, SymNum(b))
Base.:*(a::Number, b::AbstractSymbolic) = SymExpr(:*, SymNum(a), b)

# Division
Base.:/(a::AbstractSymbolic, b::AbstractSymbolic) = SymExpr(:/, a, b)
Base.:/(a::AbstractSymbolic, b::Number) = SymExpr(:/, a, SymNum(b))
Base.:/(a::Number, b::AbstractSymbolic) = SymExpr(:/, SymNum(a), b)
Base.://(a::AbstractSymbolic, b::AbstractSymbolic) = SymExpr(:/, a, b)
Base.://(a::AbstractSymbolic, b::Number) = SymExpr(:/, a, SymNum(b))
Base.://(a::Number, b::AbstractSymbolic) = SymExpr(:/, SymNum(a), b)

# Power
Base.:^(a::AbstractSymbolic, b::AbstractSymbolic) = SymExpr(:^, a, b)
Base.:^(a::AbstractSymbolic, b::Number) = SymExpr(:^, a, SymNum(b))
Base.:^(a::AbstractSymbolic, b::Integer) = SymExpr(:^, a, SymNum(b))
Base.:^(a::Number, b::AbstractSymbolic) = SymExpr(:^, SymNum(a), b)

# Square root (√ is already an alias for sqrt in Julia)
Base.sqrt(a::AbstractSymbolic) = SymExpr(:sqrt, a)

# Complex conjugate - respects :real assumption
Base.conj(a::Sym) = is_real(a) ? a : SymExpr(:conj, a)
Base.conj(a::SymNum) = SymNum(conj(a.value))
Base.conj(a::AbstractSymbolic) = SymExpr(:conj, a)

Base.adjoint(a::Sym) = conj(a)
Base.adjoint(a::SymNum) = conj(a)
Base.adjoint(a::AbstractSymbolic) = conj(a)

# Real/Imag parts - respect assumptions
Base.real(a::Sym) = is_real(a) ? a : SymExpr(:real, a)
Base.real(a::SymNum) = SymNum(real(a.value))
Base.real(a::AbstractSymbolic) = SymExpr(:real, a)

Base.imag(a::Sym) = is_real(a) ? SymNum(0) : SymExpr(:imag, a)
Base.imag(a::SymNum) = SymNum(imag(a.value))
Base.imag(a::AbstractSymbolic) = SymExpr(:imag, a)

# Absolute value
Base.abs(a::AbstractSymbolic) = SymExpr(:abs, a)
Base.abs2(a::AbstractSymbolic) = SymExpr(:abs2, a)

# Trigonometric
Base.sin(a::AbstractSymbolic) = SymExpr(:sin, a)
Base.cos(a::AbstractSymbolic) = SymExpr(:cos, a)
Base.tan(a::AbstractSymbolic) = SymExpr(:tan, a)
Base.exp(a::AbstractSymbolic) = SymExpr(:exp, a)
Base.log(a::AbstractSymbolic) = SymExpr(:log, a)

# ============== Evaluation ==============

@doc """
    substitute(expr, pairs...)
    substitute(expr, dict)

Substitute symbolic variables with values.

# Examples
```julia
n = Sym(:n)
expr = √n + 1
substitute(expr, :n => 4)  # → √4 + 1 (still symbolic with SymNum)
```
""" substitute
function substitute(s::Sym, subs::Dict{Symbol, <:Any})
    haskey(subs, s.name) ? _wrap(subs[s.name]) : s
end

substitute(s::SymNum, subs::Dict{Symbol, <:Any}) = s

function substitute(e::SymExpr, subs::Dict{Symbol, <:Any})
    new_args = map(a -> substitute(a, subs), e.args)
    SymExpr(e.op, new_args...)
end

substitute(x::Number, subs::Dict{Symbol, <:Any}) = x

# Convenience: substitute(expr, :x => 1, :y => 2)
substitute(expr, pairs::Pair{Symbol}...) = substitute(expr, Dict(pairs...))

_wrap(x::AbstractSymbolic) = x
_wrap(x::Number) = SymNum(x)

@doc """
    evaluate(expr)

Evaluate a symbolic expression to a numeric value.
All symbols must have been substituted.

# Examples
```julia
n = Sym(:n)
expr = √n + 1
result = substitute(expr, :n => 4)
evaluate(result)  # → 3.0
```
""" evaluate
evaluate(s::Sym) = error("Cannot evaluate: unsubstituted symbol $(s.name)")
evaluate(s::SymNum) = s.value
evaluate(x::Number) = x

function evaluate(e::SymExpr)
    vals = map(evaluate, e.args)
    op = e.op
    
    if op == :+
        return sum(vals)
    elseif op == :-
        return vals[1] - vals[2]
    elseif op == :neg
        return -vals[1]
    elseif op == :*
        return prod(vals)
    elseif op == :/
        return vals[1] / vals[2]
    elseif op == :^
        return vals[1] ^ vals[2]
    elseif op == :sqrt
        return sqrt(vals[1])
    elseif op == :conj
        return conj(vals[1])
    elseif op == :abs
        return abs(vals[1])
    elseif op == :abs2
        return abs2(vals[1])
    elseif op == :sin
        return sin(vals[1])
    elseif op == :cos
        return cos(vals[1])
    elseif op == :tan
        return tan(vals[1])
    elseif op == :exp
        return exp(vals[1])
    elseif op == :log
        return log(vals[1])
    else
        error("Unknown operation: $op")
    end
end

@doc """
    symbols(expr)

Get all symbolic variables in an expression.
""" symbols
symbols(s::Sym) = Set([s.name])
symbols(::SymNum) = Set{Symbol}()
symbols(e::SymExpr) = union(map(symbols, e.args)...)
symbols(::Number) = Set{Symbol}()

@doc """
    is_numeric(expr)

Check if expression contains no symbolic variables.
""" is_numeric
is_numeric(s::Sym) = false
is_numeric(::SymNum) = true
is_numeric(e::SymExpr) = all(is_numeric, e.args)
is_numeric(::Number) = true

# ============== Simplification (basic) ==============

@doc """
    simplify(expr)

Apply basic algebraic simplifications.
""" simplify
simplify(s::Sym) = s
simplify(s::SymNum) = s
simplify(x::Number) = SymNum(x)

function simplify(e::SymExpr)
    # First simplify arguments
    args = map(simplify, e.args)
    op = e.op
    
    # Identity simplifications (before numeric evaluation)
    if op == :+ && length(args) == 2
        if args[1] isa SymNum && args[1].value == 0
            return args[2]
        elseif args[2] isa SymNum && args[2].value == 0
            return args[1]
        end
    elseif op == :* && length(args) == 2
        if args[1] isa SymNum && args[1].value == 1
            return args[2]
        elseif args[2] isa SymNum && args[2].value == 1
            return args[1]
        elseif args[1] isa SymNum && args[1].value == 0
            return SymNum(0)
        elseif args[2] isa SymNum && args[2].value == 0
            return SymNum(0)
        end
    elseif op == :^ && length(args) == 2
        if args[2] isa SymNum && args[2].value == 0
            return SymNum(1)
        elseif args[2] isa SymNum && args[2].value == 1
            return args[1]
        end
    elseif op == :sqrt
        if args[1] isa SymNum && args[1].value == 0
            return SymNum(0)
        elseif args[1] isa SymNum && args[1].value == 1
            return SymNum(1)
        end
    end
    
    # If all arguments are purely numeric (SymNum only, not Sym or SymExpr), evaluate fully
    if all(a -> a isa SymNum, args)
        return SymNum(evaluate(SymExpr(op, args...)))
    end
    
    SymExpr(op, args...)
end

# ============== Comparison (for testing) ==============

# Symbolic equality (structural) - defined earlier for Sym
# Base.:(==)(a::Sym, b::Sym) already defined near struct definition
Base.:(==)(a::SymNum, b::SymNum) = a.value == b.value
Base.:(==)(a::SymNum, b::Number) = a.value == b
Base.:(==)(a::Number, b::SymNum) = a == b.value

function Base.:(==)(a::SymExpr, b::SymExpr)
    a.op == b.op && length(a.args) == length(b.args) && all(a.args .== b.args)
end

# SymExpr is never equal to a simple Sym or SymNum (structurally)
Base.:(==)(::SymExpr, ::Sym) = false
Base.:(==)(::Sym, ::SymExpr) = false
Base.:(==)(::SymExpr, ::SymNum) = false
Base.:(==)(::SymNum, ::SymExpr) = false
Base.:(==)(::Sym, ::SymNum) = false
Base.:(==)(::SymNum, ::Sym) = false

# Comparing with plain numbers
Base.:(==)(::SymExpr, ::Number) = false
Base.:(==)(::Number, ::SymExpr) = false
Base.:(==)(::Sym, ::Number) = false
Base.:(==)(::Number, ::Sym) = false

# Approximate equality (evaluate if possible)
function Base.isapprox(a::AbstractSymbolic, b::Number; kwargs...)
    is_numeric(a) && isapprox(evaluate(a), b; kwargs...)
end
function Base.isapprox(a::Number, b::AbstractSymbolic; kwargs...)
    is_numeric(b) && isapprox(a, evaluate(b); kwargs...)
end
function Base.isapprox(a::AbstractSymbolic, b::AbstractSymbolic; kwargs...)
    is_numeric(a) && is_numeric(b) && isapprox(evaluate(a), evaluate(b); kwargs...)
end

# ============== Type conversions ==============

Base.zero(::Type{<:AbstractSymbolic}) = SymNum(0)
Base.one(::Type{<:AbstractSymbolic}) = SymNum(1)
Base.zero(::AbstractSymbolic) = SymNum(0)
Base.one(::AbstractSymbolic) = SymNum(1)

# Allow symbolic in numeric contexts where needed
Base.promote_rule(::Type{<:AbstractSymbolic}, ::Type{T}) where {T<:Number} = AbstractSymbolic
Base.convert(::Type{AbstractSymbolic}, x::Number) = SymNum(x)

# ============== Kronecker Delta ==============

@doc """
    KroneckerDelta(i, j)

Represents the Kronecker delta δᵢⱼ which equals 1 if i==j and 0 otherwise.
Used for symbolic inner products where indices may be symbolic.

# Examples
```julia
k, l = Sym(:k), Sym(:l)
δ = KroneckerDelta(k, l)  # δₖₗ
evaluate(δ, :k => 2, :l => 2)  # → 1
evaluate(δ, :k => 2, :l => 3)  # → 0
```
""" KroneckerDelta
struct KroneckerDelta <: AbstractSymbolic
    i::Any  # Can be Symbol, AbstractSymbolic, or concrete value
    j::Any
end

# Display
function Base.show(io::IO, δ::KroneckerDelta)
    print(io, "δ(", δ.i, ",", δ.j, ")")
end

# Simplify KroneckerDelta - try to evaluate if indices are concrete
function simplify(δ::KroneckerDelta)
    # Simplify indices first
    i_new = δ.i isa SymExpr ? simplify(δ.i) : δ.i
    j_new = δ.j isa SymExpr ? simplify(δ.j) : δ.j
    new_δ = KroneckerDelta(i_new, j_new)
    result = _try_evaluate_delta(new_δ)
    result !== nothing ? SymNum(result) : new_δ
end

# Evaluate when both indices are concrete
function _try_evaluate_delta(δ::KroneckerDelta)
    i, j = δ.i, δ.j
    # Try to simplify symbolic expressions first
    i_simplified = i isa SymExpr ? simplify(i) : i
    j_simplified = j isa SymExpr ? simplify(j) : j
    
    # If both are concrete (not symbolic or fully numeric), we can evaluate
    i_concrete = !(i_simplified isa AbstractSymbolic) || (i_simplified isa SymNum)
    j_concrete = !(j_simplified isa AbstractSymbolic) || (j_simplified isa SymNum)
    
    if i_concrete && j_concrete
        i_val = _to_comparable_value(i_simplified)
        j_val = _to_comparable_value(j_simplified)
        return i_val == j_val ? 1 : 0
    end
    return nothing  # Cannot evaluate yet
end

# Helper to convert index to comparable value
function _to_comparable_value(x)
    if x isa AbstractSymbolic
        return evaluate(x)
    elseif x isa Symbol
        # Try to parse symbol as number (e.g., :0 → 0, :42 → 42)
        str = string(x)
        parsed = tryparse(Int, str)
        return parsed !== nothing ? parsed : x
    else
        return x
    end
end

# Check if delta is numerically evaluable
is_numeric(δ::KroneckerDelta) = _try_evaluate_delta(δ) !== nothing

# Evaluate delta
function evaluate(δ::KroneckerDelta)
    result = _try_evaluate_delta(δ)
    result !== nothing && return result
    error("Cannot evaluate KroneckerDelta with symbolic indices: δ($(δ.i), $(δ.j))")
end

# Substitute in delta (Pair version)
function substitute(δ::KroneckerDelta, replacements::Pair{Symbol}...)
    i_new = δ.i isa AbstractSymbolic ? substitute(δ.i, replacements...) : δ.i
    j_new = δ.j isa AbstractSymbolic ? substitute(δ.j, replacements...) : δ.j
    
    # Simplify after substitution
    i_simplified = i_new isa SymExpr ? simplify(i_new) : i_new
    j_simplified = j_new isa SymExpr ? simplify(j_new) : j_new
    
    # Try to evaluate delta
    new_δ = KroneckerDelta(i_simplified, j_simplified)
    result = _try_evaluate_delta(new_δ)
    result !== nothing ? result : new_δ
end

# Substitute in delta (Dict version - called from SymExpr substitute)
function substitute(δ::KroneckerDelta, subs::Dict{Symbol, <:Any})
    i_new = δ.i isa AbstractSymbolic ? substitute(δ.i, subs) : δ.i
    j_new = δ.j isa AbstractSymbolic ? substitute(δ.j, subs) : δ.j
    
    # Simplify after substitution
    i_simplified = i_new isa SymExpr ? simplify(i_new) : i_new
    j_simplified = j_new isa SymExpr ? simplify(j_new) : j_new
    
    # Try to evaluate delta
    new_δ = KroneckerDelta(i_simplified, j_simplified)
    result = _try_evaluate_delta(new_δ)
    result !== nothing ? result : new_δ
end

# Arithmetic with KroneckerDelta
Base.:*(a::Number, δ::KroneckerDelta) = iszero(a) ? 0 : SymExpr(:*, SymNum(a), δ)
Base.:*(δ::KroneckerDelta, a::Number) = a * δ
Base.:*(a::AbstractSymbolic, δ::KroneckerDelta) = SymExpr(:*, a, δ)
Base.:*(δ::KroneckerDelta, a::AbstractSymbolic) = SymExpr(:*, δ, a)
Base.:*(δ1::KroneckerDelta, δ2::KroneckerDelta) = SymExpr(:*, δ1, δ2)

Base.:+(δ::KroneckerDelta, a::Number) = SymExpr(:+, δ, SymNum(a))
Base.:+(a::Number, δ::KroneckerDelta) = SymExpr(:+, SymNum(a), δ)
Base.:+(δ::KroneckerDelta, a::AbstractSymbolic) = SymExpr(:+, δ, a)
Base.:+(a::AbstractSymbolic, δ::KroneckerDelta) = SymExpr(:+, a, δ)

# Comparisons for delta
Base.iszero(δ::KroneckerDelta) = _try_evaluate_delta(δ) == 0
Base.isone(δ::KroneckerDelta) = _try_evaluate_delta(δ) == 1
Base.:(==)(δ1::KroneckerDelta, δ2::KroneckerDelta) = δ1.i == δ2.i && δ1.j == δ2.j
