
using LinearAlgebra
using SparseArrays
using SuiteSparse

@testset "Interface with LinearAlgebra" begin
    n = 512
    M = sprandn(n, n, .2)
    b = randn(n)
    # TODO: add diagonal term to avoid nonzero F in KLU
    A0 = M + 1e-4I
    A1 = M +    2I
    #= End of MIT LICENSE =#

    sol1 = A0  \ b
    sol2 = A0' \ b
    sol3 = A1  \ b

    A0_d = CuSparseMatrixCSR(A0)
    A1_d = CuSparseMatrixCSR(A1)

    @testset "Normal mode ($sym)" for sym in [:KLU]
        x_d = CUDA.zeros(Float64, n)
        b_d = CuVector(b)
        rf = CUSOLVERRF.RFLU(A0; symbolic=sym)

        ldiv!(x_d, rf, b_d)
        @test Vector(x_d) ≈ sol1

        ldiv!(x_d, rf', b_d)
        @test Vector(x_d) ≈ sol2

        # Refactorization
        lu!(rf, A1_d)

        ldiv!(x_d, rf, b_d)
        @test Vector(x_d) ≈ sol3
    end

    @testset "Batch mode ($sym)" for sym in [:KLU]
        nrhs = 64
        B = rand(n, nrhs)
        X0_sol = A0  \ B
        X1_sol = A0' \ B
        X2_sol = A1  \ B

        X_d = CUDA.zeros(Float64, n, nrhs)
        B_d = B |> CuMatrix

        rf = CUSOLVERRF.RFLU(A0; symbolic=sym, nrhs=nrhs)

        ldiv!(X_d, rf, B_d)
        @test Matrix(X_d) ≈ X0_sol

        ldiv!(X_d, rf', B_d)
        @test Matrix(X_d) ≈ X1_sol

        # Refactorization
        lu!(rf, A1_d)

        ldiv!(X_d, rf, B_d)
        @test Matrix(X_d) ≈ X2_sol
    end
end

