using TestItemRunner
@run_package_tests

@testitem "Hilbert Spaces" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    @test H == HilbertSpace(:H, 2)
    @test H != HilbertSpace(:H, 3)
    @test H != HilbertSpace(:G, 2)
    
    F = FockSpace(:F)
    @test F == HilbertSpace(:F, nothing)
    
    # Composite spaces
    H1 = HilbertSpace(:A, 2)
    H2 = HilbertSpace(:B, 2)
    H12 = H1 ⊗ H2
    @test H12 == CompositeSpace(H1, H2)
end

@testitem "Basis" begin
    using QSymbolic
    
    H = HilbertSpace(:spin, 2)
    Zb = Basis(H, :z)
    Xb = Basis(H, :x)
    
    @test Zb == Basis(H, :z)
    @test Zb != Xb
    @test space(Zb) == typeof(H)
    @test basisname(Zb) == :z
end

@testitem "BasisKet and BasisBra" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    
    # Default basis (backward compatible)
    ψ = BasisKet(H, :ψ)
    ϕ = BasisKet(H, :ϕ)
    
    @test ψ == BasisKet(H, :ψ)
    @test ψ != ϕ
    
    # Adjoint
    @test ψ' isa BasisBra
    @test (ψ')' == ψ
    
    # Inner products (same basis → orthonormal)
    @test ψ' * ψ == 1
    @test ψ' * ϕ == 0
    
    # Explicit basis
    Zb = Basis(H, :z)
    up = BasisKet(Zb, :↑)
    down = BasisKet(Zb, :↓)
    @test up' * up == 1
    @test up' * down == 0
end

@testitem "Cross-basis inner products" begin
    using QSymbolic
    
    H = HilbertSpace(:spin, 2)
    Zb = Basis(H, :z)
    Xb = Basis(H, :x)
    
    up_z = BasisKet(Zb, :↑)
    up_x = BasisKet(Xb, :↑)
    
    # No transform defined → symbolic
    result = up_z' * up_x
    @test result isa InnerProduct
    
    # Define transform: |↑ₓ⟩ = (|↑ᵤ⟩ + |↓ᵤ⟩)/√2
    down_z = BasisKet(Zb, :↓)
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
    ψ = BasisKet(H, :ψ)
    
    wket = 2 * ψ
    @test wket isa weightedKet
    @test wket.weight == 2
    
    wket2 = 3 * wket
    @test wket2.weight == 6
    
    # Division
    @test (ψ / 2).weight == 0.5
    @test (ψ // 2).weight == 1//2
    
    # Bras
    wbra = 2 * ψ'
    @test wbra isa weightedBra
end

@testitem "Addition" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    ψ = BasisKet(H, :ψ)
    ϕ = BasisKet(H, :ϕ)
    
    sum_ket = ψ + ϕ
    @test sum_ket isa sumKet
    @test length(sum_ket.kets) == 2
    
    # Subtraction
    diff_ket = ψ - ϕ
    @test diff_ket isa sumKet
    @test diff_ket.weights == [1, -1]
    
    # Same ket addition
    double_ψ = ψ + ψ
    @test double_ψ isa weightedKet
    @test double_ψ.weight == 2
end

@testitem "Inner products with sums" begin
    using QSymbolic
    
    H = HilbertSpace(:H, 2)
    ψ = BasisKet(H, :ψ)
    ϕ = BasisKet(H, :ϕ)
    
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
    H2 = HilbertSpace(:B, 2)
    
    ψ = BasisKet(H1, :ψ)
    ϕ = BasisKet(H2, :ϕ)
    
    # Tensor product
    product = ψ ⊗ ϕ
    @test product isa ProductKet
    
    # Adjoint
    @test product' isa ProductBra
    @test (product')' == product
    
    # Inner product (same composite basis)
    @test product' * product == 1
    
    # Orthogonal states
    χ = BasisKet(H2, :χ)
    product2 = ψ ⊗ χ
    @test product' * product2 == 0
end

@testitem "Composite basis transforms" begin
    using QSymbolic
    
    H1 = HilbertSpace(:A, 2)
    H2 = HilbertSpace(:B, 2)
    
    # Bases for subsystem A
    Za = Basis(H1, :z)
    Xa = Basis(H1, :x)
    up_a = BasisKet(Za, :↑)
    down_a = BasisKet(Za, :↓)
    
    # Bases for subsystem B  
    Zb = Basis(H2, :z)
    Xb = Basis(H2, :x)
    up_b = BasisKet(Zb, :↑)
    down_b = BasisKet(Zb, :↓)
    
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
    up_x_a = BasisKet(Xa, :↑)
    state = up_x_a ⊗ up_b  # |↑ₓ⟩_A ⊗ |↑ᵤ⟩_B
    
    # Transform to z⊗z basis uses factorized transform
    transformed = transform(state, typeof(ZaZb))
    @test transformed isa SumProductKet
    @test length(transformed.kets) == 2  # (|↑⟩+|↓⟩)/√2 ⊗ |↑⟩
    
    # Verify inner product: ⟨↑↑|transformed⟩ = 1/√2
    zz_up = up_a ⊗ up_b
    @test zz_up' * transformed ≈ 1/√2
    
    clear_transforms!()
end

@testitem "Symbolic scalars" begin
    using QSymbolic
    
    # Create symbolic variables
    n = Sym(:n)
    m = Sym(:m)
    
    @test n isa AbstractSymbolic
    @test n == Sym(:n)
    @test n != m
    
    # Arithmetic builds expressions
    @test √n isa SymExpr
    @test (n + 1) isa SymExpr
    @test (n * m) isa SymExpr
    @test (n / 2) isa SymExpr
    @test (n^2) isa SymExpr
    
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

@testitem "Symbolic substitution and evaluation" begin
    using QSymbolic
    
    n = Sym(:n)
    
    # Substitute
    expr = √n + 1
    result = substitute(expr, :n => 4)
    @test is_numeric(result)
    @test evaluate(result) ≈ 3.0
    
    # Multiple substitutions
    m = Sym(:m)
    expr2 = n * m
    result2 = substitute(expr2, :n => 2, :m => 3)
    @test evaluate(result2) == 6
    
    # Partial substitution
    partial = substitute(expr2, :n => 2)
    @test !is_numeric(partial)
    @test :m in symbols(partial)
end

@testitem "Symbolic simplification" begin
    using QSymbolic
    
    n = Sym(:n)
    
    # Identity simplifications
    @test simplify(n * 1) == n
    @test simplify(n * 0) == SymNum(0)
    @test simplify(n + 0) == n
    @test simplify(n^1) == n
    @test simplify(n^0) == SymNum(1)
    
    # Numeric folding
    @test simplify(SymNum(2) + SymNum(3)) == SymNum(5)
    @test simplify(SymNum(2) * SymNum(3)) == SymNum(6)
end
