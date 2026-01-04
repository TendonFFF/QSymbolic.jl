# Custom contraction rules for basis inner products
# Users can define specific rules for how indices contract in different bases

# Global dictionary to store custom contraction rules
# Key: Type{<:AbstractBasis}
# Value: Function that takes (bra_index, ket_index) and returns result
const CONTRACTION_RULES = Dict{Type{<:AbstractBasis}, Function}()

@doc """
    define_contraction_rule!(basis::Type{<:AbstractBasis}, rule::Function)
    define_contraction_rule!(basis_type, rule) do bra_idx, ket_idx
        # Custom contraction logic
        return result
    end

Define a custom contraction rule for computing inner products in a specific basis.

The rule function should take two arguments:
- `bra_idx`: The index of the bra (can be single value or tuple)
- `ket_idx`: The index of the ket (can be single value or tuple)

And return:
- A number (scalar result)
- A symbolic expression (e.g., KroneckerDelta)
- A Dirac delta expression

# Examples
```julia
# Define a contraction rule for Fock basis that returns overlap integrals
F = FockSpace(:F)
Fb = Basis(F, :n)

define_contraction_rule!(typeof(Fb)) do i, j
    # For Fock states, orthonormal: ⟨i|j⟩ = δᵢⱼ
    if i isa AbstractSymbolic || j isa AbstractSymbolic
        return isequal(i, j) ? 1 : KroneckerDelta(i, j)
    else
        return i == j ? 1 : 0
    end
end

# Define a non-orthonormal basis with overlap
H = HilbertSpace(:H, 2)
Cb = Basis(H, :coherent)

@syms α β
define_contraction_rule!(typeof(Cb)) do i, j
    # Coherent states: ⟨α|β⟩ = exp(-|α-β|²/2)
    if i isa AbstractSymbolic && j isa AbstractSymbolic
        # Return symbolic expression
        exp(-abs(i - j)^2 / 2)
    else
        # Evaluate numerically
        exp(-abs(i - j)^2 / 2)
    end
end
```

See also: [`has_contraction_rule`](@ref), [`clear_contraction_rules!`](@ref)
"""
function define_contraction_rule!(basis::Type{<:AbstractBasis}, rule::Function)
    CONTRACTION_RULES[basis] = rule
    return nothing
end

# Convenience version with do-block syntax
function define_contraction_rule!(f::Function, basis::Type{<:AbstractBasis})
    define_contraction_rule!(basis, f)
end

@doc """
    has_contraction_rule(basis::Type{<:AbstractBasis}) -> Bool

Check if a custom contraction rule has been defined for the given basis.

# Examples
```julia
H = HilbertSpace(:H, 2)
Zb = Basis(H, :z)

has_contraction_rule(typeof(Zb))  # false (uses default)

define_contraction_rule!(typeof(Zb)) do i, j
    i == j ? 1 : 0
end

has_contraction_rule(typeof(Zb))  # true
```

See also: [`define_contraction_rule!`](@ref)
"""
function has_contraction_rule(basis::Type{<:AbstractBasis})
    return haskey(CONTRACTION_RULES, basis)
end

@doc """
    get_contraction_rule(basis::Type{<:AbstractBasis}) -> Function

Get the custom contraction rule for a basis, or return nothing if not defined.

See also: [`define_contraction_rule!`](@ref), [`has_contraction_rule`](@ref)
"""
function get_contraction_rule(basis::Type{<:AbstractBasis})
    return get(CONTRACTION_RULES, basis, nothing)
end

@doc """
    clear_contraction_rules!()

Clear all custom contraction rules, reverting to default index-based checking.

See also: [`define_contraction_rule!`](@ref)
"""
function clear_contraction_rules!()
    empty!(CONTRACTION_RULES)
    return nothing
end

@doc """
    apply_contraction_rule(basis::Type{<:AbstractBasis}, bra_idx, ket_idx)

Apply the contraction rule for a basis (custom if defined, otherwise default).

The default rule checks if indices are equal:
- If both concrete: return 1 if equal, 0 if not
- If any symbolic: return 1 if symbolically equal, KroneckerDelta otherwise
- For tuples: check component-wise

# Examples
```julia
H = HilbertSpace(:H, 2)
Zb = Basis(H, :z)

# Default rule (orthonormal)
apply_contraction_rule(typeof(Zb), :↑, :↑)  # 1
apply_contraction_rule(typeof(Zb), :↑, :↓)  # 0

# With symbolic
n = Sym(:n)
apply_contraction_rule(typeof(Zb), n, n)  # 1
apply_contraction_rule(typeof(Zb), n, :↑)  # KroneckerDelta(n, :↑)
```

See also: [`define_contraction_rule!`](@ref)
"""
function apply_contraction_rule(basis::Type{<:AbstractBasis}, bra_idx, ket_idx)
    # Check if custom rule exists
    if has_contraction_rule(basis)
        rule = get_contraction_rule(basis)
        return rule(bra_idx, ket_idx)
    end
    
    # Default rule: orthonormal basis with index checking
    return _default_contraction_rule(bra_idx, ket_idx)
end

# Default contraction rule for orthonormal bases
function _default_contraction_rule(bra_idx, ket_idx)
    # Handle multi-index case: must check all components
    if bra_idx isa Tuple && ket_idx isa Tuple
        # Both are tuples - check if same length
        if length(bra_idx) != length(ket_idx)
            # Different number of indices → orthogonal
            return 0
        end
        
        # Check if any index is symbolic
        has_symbolic = any(idx isa AbstractSymbolic for idx in bra_idx) || 
                      any(idx isa AbstractSymbolic for idx in ket_idx)
        
        if has_symbolic
            # If any component is symbolic, need Kronecker delta
            # Check if all components are symbolically equal
            if all(isequal(bra_idx[k], ket_idx[k]) for k in 1:length(bra_idx))
                return 1
            else
                # Return product of Kronecker deltas for each component
                result = 1
                for k in 1:length(bra_idx)
                    if bra_idx[k] isa AbstractSymbolic || ket_idx[k] isa AbstractSymbolic
                        if isequal(bra_idx[k], ket_idx[k])
                            # Symbolically equal component contributes 1
                            continue
                        else
                            result = simplify(result * KroneckerDelta(bra_idx[k], ket_idx[k]))
                        end
                    else
                        # Concrete component
                        if bra_idx[k] != ket_idx[k]
                            return 0  # Orthogonal
                        end
                    end
                end
                return result
            end
        else
            # All concrete: evaluate directly
            return all(bra_idx[k] == ket_idx[k] for k in 1:length(bra_idx)) ? 1 : 0
        end
    elseif bra_idx isa Tuple || ket_idx isa Tuple
        # One is tuple, one is not → incompatible indices
        return 0
    else
        # Single indices (not tuples)
        if bra_idx isa AbstractSymbolic || ket_idx isa AbstractSymbolic
            # Try to simplify if indices are symbolically equal
            if isequal(bra_idx, ket_idx)
                return 1
            else
                return KroneckerDelta(bra_idx, ket_idx)
            end
        else
            # Both concrete: evaluate directly
            return bra_idx == ket_idx ? 1 : 0
        end
    end
end
