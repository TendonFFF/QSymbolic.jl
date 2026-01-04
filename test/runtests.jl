using TestItemRunner
@run_package_tests

@testitem "Hilbert Spaces" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    Hb = Basis(H, :default)
    @test H == HilbertSpace(:H, 2)
    @test H != HilbertSpace(:H, 3)
    @test H != HilbertSpace(:G, 2)
    
    F = FockSpace(:F)
    @test F == HilbertSpace(:F, nothing)
    
    # Composite spaces
    H1 = HilbertSpace(:A, 2)
    H1b = Basis(H1, :default)
    H2 = HilbertSpace(:B, 2)
    H2b = Basis(H2, :default)
    H12 = H1 ⊗ H2
    @test H12 == CompositeSpace(H1, H2)
end

@testitem "Basis" begin
    using QSymbolic
    
    H = HilbertSpace(:spin, 2)
    Hb = Basis(H, :default)
    Zb = Basis(H, :z)
    Xb = Basis(H, :x)
    
    @test Zb == Basis(H, :z)
    @test Zb != Xb
    @test space(Zb) == typeof(H)
    @test basisname(Zb) == :z
end

@testitem "Ket and Bra" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    Hb = Basis(H, :default)
    
    # Default basis (backward compatible)
    ψ = Ket(Hb, :ψ)
    ϕ = Ket(Hb, :ϕ)
    
    @test ψ == Ket(Hb, :ψ)
    @test ψ != ϕ
    
    # Adjoint
    @test ψ' isa Bra
    @test (ψ')' == ψ
    
    # Inner products (same basis → orthonormal)
    @test ψ' * ψ == 1
    @test ψ' * ϕ == 0
    
    # Explicit basis
    Zb = Basis(H, :z)
    up = Ket(Zb, :↑)
    down = Ket(Zb, :↓)
    @test up' * up == 1
    @test up' * down == 0
end

@testitem "Cross-basis inner products" begin
    using QSymbolic
    
    H = HilbertSpace(:spin, 2)
    Hb = Basis(H, :default)
    Zb = Basis(H, :z)
    Xb = Basis(H, :x)
    
    up_z = Ket(Zb, :↑)
    up_x = Ket(Xb, :↑)
    
    # No transform defined → symbolic
    result = up_z' * up_x
    @test result isa InnerProduct
    
    # Define transform: |↑ₓ⟩ = (|↑ᵤ⟩ + |↓ᵤ⟩)/√2
    down_z = Ket(Zb, :↓)
    define_transform!(Xb, Zb) do idx
        if idx == :↑
            (up_z + down_z) / √2
        else
            (up_z - down_z) / √2
        end
    end
    
    # Now cross-basis works
    @test up_z' * up_x ≈ 1/√2
    @test down_z' * up_x ≈ 1/√2
    
    clear_transforms!()
end

@testitem "Scalar multiplication" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    Hb = Basis(H, :default)
    ψ = Ket(Hb, :ψ)
    
    wket = 2 * ψ
    @test wket isa WeightedKet
    @test wket.weight == 2
    
    wket2 = 3 * wket
    @test wket2.weight == 6
    
    # Division
    @test (ψ / 2).weight == 0.5
    @test (ψ // 2).weight == 1//2
    
    # Bras
    wbra = 2 * ψ'
    @test wbra isa WeightedBra
end

@testitem "Addition" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    Hb = Basis(H, :default)
    ψ = Ket(Hb, :ψ)
    ϕ = Ket(Hb, :ϕ)
    
    sum_ket = ψ + ϕ
    @test sum_ket isa SumKet
    @test length(sum_ket.kets) == 2
    
    # Subtraction
    diff_ket = ψ - ϕ
    @test diff_ket isa SumKet
    @test diff_ket.weights == [1, -1]
    
    # Same ket addition
    double_ψ = ψ + ψ
    @test double_ψ isa WeightedKet
    @test double_ψ.weight == 2
end

@testitem "Inner products with sums" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    Hb = Basis(H, :default)
    ψ = Ket(Hb, :ψ)
    ϕ = Ket(Hb, :ϕ)
    
    # ⟨ψ|ψ⟩ + ⟨ψ|ϕ⟩ = 1 + 0 = 1
    @test (ψ + ϕ)' * ψ == 1
    
    # (⟨ψ| + ⟨ϕ|)(|ψ⟩ + |ϕ⟩) = 1 + 0 + 0 + 1 = 2
    @test (ψ + ϕ)' * (ψ + ϕ) == 2
end

@testitem "Fock space" begin
    using QSymbolic
    
    F = FockSpace(:F)
    n0 = FockKet(F, 0)
    n1 = FockKet(F, 1)
    
    @test n0' * n0 == 1
    @test n0' * n1 == 0
end

@testitem "Product states" begin
    using QSymbolic
    
    H1 = HilbertSpace(:A, 2)
    H1b = Basis(H1, :default)
    H2 = HilbertSpace(:B, 2)
    H2b = Basis(H2, :default)
    
    ψ = Ket(H1b, :ψ)
    ϕ = Ket(H2b, :ϕ)
    
    # Tensor product
    product = ψ ⊗ ϕ
    @test product isa ProductKet
    
    # Adjoint
    @test product' isa ProductBra
    @test (product')' == product
    
    # Inner product (same composite basis)
    @test product' * product == 1
    
    # Orthogonal states
    χ = Ket(H2b, :χ)
    product2 = ψ ⊗ χ
    @test product' * product2 == 0
end

@testitem "Composite basis transforms" begin
    using QSymbolic
    
    H1 = HilbertSpace(:A, 2)
    H1b = Basis(H1, :default)
    H2 = HilbertSpace(:B, 2)
    H2b = Basis(H2, :default)
    
    # Bases for subsystem A
    Za = Basis(H1, :z)
    Xa = Basis(H1, :x)
    up_a = Ket(Za, :↑)
    down_a = Ket(Za, :↓)
    
    # Bases for subsystem B  
    Zb = Basis(H2, :z)
    Xb = Basis(H2, :x)
    up_b = Ket(Zb, :↑)
    down_b = Ket(Zb, :↓)
    
    # Define transforms for each subsystem
    define_transform!(Xa, Za) do idx
        idx == :↑ ? (up_a + down_a) / √2 : (up_a - down_a) / √2
    end
    define_transform!(Xb, Zb) do idx
        idx == :↑ ? (up_b + down_b) / √2 : (up_b - down_b) / √2
    end
    
    # Composite bases
    ZaZb = Za ⊗ Zb
    XaZb = Xa ⊗ Zb
    
    @test ZaZb isa CompositeBasis
    @test has_transform(typeof(XaZb), typeof(ZaZb))
    
    # Product state in x⊗z basis
    up_x_a = Ket(Xa, :↑)
    state = up_x_a ⊗ up_b  # |↑ₓ⟩_A ⊗ |↑ᵤ⟩_B
    
    # Transform to z⊗z basis uses factorized transform
    transformed = transform(state, typeof(ZaZb))
    @test transformed isa SumKet
    @test length(transformed.kets) == 2  # (|↑⟩+|↓⟩)/√2 ⊗ |↑⟩
    
    # Verify inner product: ⟨↑↑|transformed⟩ = 1/√2
    zz_up = up_a ⊗ up_b
    @test zz_up' * transformed ≈ 1/√2
    
    clear_transforms!()
end

@testitem "Symbolic scalars" begin
    using QSymbolic
    import Symbolics
    
    # Create symbolic variables
    n = Sym(:n)
    m = Sym(:m)
    
    @test n isa AbstractSymbolic
    @test isequal(n, Sym(:n))  # Use isequal for structural equality
    @test !isequal(n, m)
    
    # Arithmetic builds expressions (Symbolics.Num)
    @test √n isa Symbolics.Num
    @test (n + 1) isa Symbolics.Num
    @test (n * m) isa Symbolics.Num
    @test (n / 2) isa Symbolics.Num
    @test (n^2) isa Symbolics.Num
    
    # Symbol extraction
    expr = n^2 + 2*n*m + m^2
    @test :n in symbols(expr)
    @test :m in symbols(expr)
    @test length(symbols(expr)) == 2
    
    # Numeric check
    @test !is_numeric(n)
    @test is_numeric(SymNum(5))
    @test !is_numeric(n + 1)
    @test is_numeric(SymNum(2) + SymNum(3))
end

@testitem "Symbolic assumptions" begin
    using QSymbolic
    import Symbolics
    
    # Generic symbol - Symbolics.jl defaults to Real
    z = Sym(:z)
    @test is_real(z)  # Symbolics.jl symbols are Real by default
    @test !is_positive(z)
    @test !is_integer(z)
    @test :real in assumptions(z)  # Default is Real in Symbolics.jl
    
    # Real symbol using positional args
    r = Sym(:r, :real)
    @test is_real(r)
    @test :real in assumptions(r)
    
    # Positive - Note: Symbolics.jl doesn't track positive as a type constraint
    # We store it as Real (Symbolics.jl limitation)
    p = Sym(:p, :positive)
    @test is_real(p)  # positive implies real
    @test :real in assumptions(p)
    # is_positive and is_nonnegative return false for symbolic (Symbolics.jl limitation)
    
    # Integer implies real
    k = Sym(:k, :integer)
    @test is_integer(k)
    @test is_real(k)
    
    # Using keyword syntax
    m = Sym(:m, integer=true, positive=true)
    @test is_integer(m)
    @test is_real(m)
    
    # Adjoint/conjugate for real symbols
    @test conj(z) isa Symbolics.Num
    @test isequal(conj(r), r)  # real: conj is identity
    @test z' isa Symbolics.Num
    @test isequal(r', r)       # real: adjoint is identity
    
    # Real part for real symbols
    @test isequal(real(r), r)  # real: real() is identity
    
    # Imag part for real is zero
    @test isequal(imag(r), 0)
    
    # Equality uses isequal for structural comparison
    @test isequal(Sym(:x), Sym(:x))
    @test !isequal(Sym(:x), Sym(:y))
end

@testitem "Symbolic substitution and evaluation" begin
    using QSymbolic
    import Symbolics
    
    n = Sym(:n)
    
    # Substitute - use linear expression for simpler evaluation
    expr = n + 1
    result = substitute(expr, :n => 2)
    # The substituted expression should be a Symbolics.Num with numeric value
    @test Float64(Symbolics.value(result)) ≈ 3.0
    
    # Multiple substitutions
    m = Sym(:m)
    expr2 = n * m
    result2 = substitute(expr2, :n => 2, :m => 3)
    @test Float64(Symbolics.value(result2)) == 6.0
    
    # Partial substitution
    partial = substitute(expr2, :n => 2)
    @test !is_numeric(partial)
    @test :m in symbols(partial)
end

@testitem "Symbolic simplification" begin
    using QSymbolic
    import Symbolics
    
    n = Sym(:n)
    
    # Identity simplifications - Symbolics.jl auto-simplifies
    @test isequal(simplify(n * 1), n)
    @test isequal(simplify(n * 0), 0)
    @test isequal(simplify(n + 0), n)
    @test isequal(simplify(n^1), n)
    @test isequal(simplify(n^0), 1)
    
    # Numeric folding
    @test simplify(SymNum(2) + SymNum(3)) == 5
    @test simplify(SymNum(2) * SymNum(3)) == 6
end

@testitem "Outer product operators" begin
    using QSymbolic
    
    H = HilbertSpace(:spin, 2)
    Hb = Basis(H, :default)
    Zb = Basis(H, :z)
    up = Ket(Zb, :↑)
    down = Ket(Zb, :↓)
    
    # Create operator via outer product
    P_up = up * up'
    @test P_up isa Operator
    
    # Projector action
    result = P_up * up
    @test result isa WeightedKet
    
    @test P_up * down == 0
    
    # Ladder operator
    σ_plus = up * down'
    @test (σ_plus * down) isa WeightedKet
    @test σ_plus * up == 0
    
    # Adjoint
    σ_minus = σ_plus'
    @test (σ_minus * up) isa WeightedKet
    @test σ_minus * down == 0
end

@testitem "Operator algebra" begin
    using QSymbolic
    
    H = HilbertSpace(:spin, 2)
    Hb = Basis(H, :default)
    Zb = Basis(H, :z)
    up = Ket(Zb, :↑)
    down = Ket(Zb, :↓)
    
    P_up = up * up'
    P_down = down * down'
    
    # Sum of operators (identity)
    I = P_up + P_down
    @test I isa SumOperator
    @test (I * up) isa WeightedKet
    @test (I * down) isa WeightedKet
    
    # Pauli Z
    σz = P_up - P_down
    @test σz isa SumOperator
    
    # Apply to superposition
    ψ = (up + down) / √2
    result = σz * ψ
    @test result isa SumKet
    
    # Operator product
    σ_plus = up * down'
    σ_minus = down * up'
    product = σ_plus * σ_minus
    @test product isa OperatorProduct
    
    # (σ+σ-)|↑⟩ = |↑⟩⟨↓|↓⟩⟨↑|↑⟩ = |↑⟩
    @test (product * up) isa WeightedKet
    @test product * down == 0
end

@testitem "Function operators" begin
    using QSymbolic
    
    F = FockSpace(:mode)
    Fb = Basis(F, :n)
    
    # Simple number operator: N|n⟩ = n|n⟩
    N = FunctionOperator(:N, Fb) do ket
        n = parse(Int, string(ket.index))
        n * ket
    end
    
    n0 = Ket(Fb, 0)
    n3 = Ket(Fb, 3)
    
    # N|0⟩ = 0|0⟩ (zero weight, not literal 0)
    result0 = N * n0
    @test result0 isa WeightedKet
    @test result0.weight == 0
    
    # N|3⟩ = 3|3⟩
    result3 = N * n3
    @test result3 isa WeightedKet
    @test result3.weight == 3
end

@testitem "FunctionOperator with adjoint action" begin
    using QSymbolic
    
    F = FockSpace(:mode)
    Fb = Basis(F, :n)
    
    # Create ladder operators with explicit adjoint
    a, adag = create_ladder_operators(Fb)
    
    n0 = Ket(Fb, 0)
    n1 = Ket(Fb, 1)
    n2 = Ket(Fb, 2)
    
    # Annihilation: a|n⟩ = √n |n-1⟩
    @test a * n0 == 0
    result = a * n1
    @test result isa WeightedKet
    @test result.weight ≈ 1.0
    @test result.Ket.index == Symbol("0")
    
    # Creation: a†|n⟩ = √(n+1) |n+1⟩
    result = adag * n0
    @test result isa WeightedKet
    @test result.weight ≈ 1.0
    @test result.Ket.index == Symbol("1")
    
    # Number operator: N = a†a
    N = adag * a
    @test N isa OperatorProduct
    
    result = N * n2
    @test result isa WeightedKet
    @test result.weight ≈ 2.0
    @test result.Ket.index == Symbol("2")
end

@testitem "FunctionOperator arithmetic" begin
    using QSymbolic
    
    F = FockSpace(:mode)
    Fb = Basis(F, :n)
    
    a, adag = create_ladder_operators(Fb)
    n0 = Ket(Fb, 0)
    n1 = Ket(Fb, 1)
    
    # Scalar multiplication
    scaled = 2 * a
    @test scaled isa ScaledOperator
    result = scaled * n1
    @test result isa WeightedKet
    @test result.weight ≈ 2.0
    
    # Addition
    sum_op = a + adag
    @test sum_op isa SumOperator
    
    # Subtraction
    diff_op = a - adag
    @test diff_op isa SumOperator
    
    # Product with other operators
    prod = a * adag
    @test prod isa OperatorProduct
end

@testitem "Symbolic indices for kets and bras" begin
    using QSymbolic
    
    F = FockSpace(:fock)
    Fb = Basis(F, :n)
    
    # Create symbolic index
    n = Sym(:n)
    
    # Ket with symbolic index
    ket_n = Ket(Fb, n)
    @test ket_n isa Ket
    @test ket_n.index isa AbstractSymbolic
    @test isequal(ket_n.index, n)
    
    # Bra with symbolic index
    bra_n = Bra(Fb, n)
    @test bra_n isa Bra
    @test bra_n.index isa AbstractSymbolic
    
    # Adjoint of ket with symbolic index
    ket_n_adj = ket_n'
    @test ket_n_adj isa Bra
    @test isequal(ket_n_adj.index, n)
    
    # Weighted ket with symbolic index
    weighted = 2 * ket_n
    @test weighted isa WeightedKet
    @test isequal(weighted.Ket.index, n)
    
    # Symbolic weights
    α = Sym(:α)
    weighted_sym = α * Ket(Fb, :0)
    @test weighted_sym isa WeightedKet
    @test weighted_sym.weight isa AbstractSymbolic
    
    # Display (should show symbolic)
    io = IOBuffer()
    show(io, ket_n)
    @test occursin("n", String(take!(io)))
end

@testitem "Symbolic arithmetic with quantum states" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    Hb = Basis(H, :default)
    ψ = Ket(Hb, :ψ)
    ϕ = Ket(Hb, :ϕ)
    
    # Symbolic coefficients in sum
    α = Sym(:α)
    β = Sym(:β)
    
    state = α * ψ + β * ϕ
    @test state isa SumKet
    @test :α in symbols(state.weights[1])
    @test :β in symbols(state.weights[2])
    
    # Substitution workflow
    concrete = substitute.(state.weights, Ref(Dict(:α => 1/√2, :β => 1/√2)))
    @test all(is_numeric, concrete)
end
@testitem "Multi-index for states" begin
    using QSymbolic
    
    # Setup composite space
    H1 = HilbertSpace(:A, 2)
    H1b = Basis(H1, :default)
    H2 = HilbertSpace(:B, 2)
    H2b = Basis(H2, :default)
    composite = H1 ⊗ H2
    
    # Multi-index with symbols
    ket_sym = Ket(composite, (:↑, :↓))
    @test ket_sym isa Ket
    @test ket_sym.index == (:↑, :↓)
    @test ket_sym.index isa Tuple
    
    # Multi-index with integers (converted to symbols)
    ket_int = Ket(composite, (0, 1))
    @test ket_int.index == (Symbol("0"), Symbol("1"))
    
    # Multi-index with symbolic
    n = Sym(:n)
    m = Sym(:m)
    ket_symbolic = Ket(composite, (n, m))
    @test ket_symbolic.index == (n, m)
    @test ket_symbolic.index[1] isa AbstractSymbolic
    @test ket_symbolic.index[2] isa AbstractSymbolic
    
    # Mixed index (symbol + symbolic)
    ket_mixed = Ket(composite, (:↑, n))
    @test ket_mixed.index == (:↑, n)
    
    # Bra multi-index
    bra_sym = Bra(composite, (:↑, :↓))
    @test bra_sym.index == (:↑, :↓)
    
    # Adjoint preserves multi-index
    bra_adj = ket_sym'
    @test bra_adj isa Bra
    @test bra_adj.index == (:↑, :↓)
    
    # Display format (comma-separated)
    io = IOBuffer()
    show(io, ket_sym)
    @test String(take!(io)) == "|↑,↓⟩"
    
    # CompositeBasis
    B1 = Basis(H1, :z)
    B2 = Basis(H2, :z)
    CB = B1 ⊗ B2
    ket_cb = Ket(CB, (n, m))
    @test ket_cb.index == (n, m)
end