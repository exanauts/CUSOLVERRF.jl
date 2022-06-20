
using LinearAlgebra
using SuiteSparse

@testset "Interface with LinearAlgebra" begin
    #=
        Code taken from https://github.com/JuliaSparse/KLU.jl/blob/main/test/runtests.jl#L9-L15
        Subject to

        MIT License

        Copyright (c) 2021 Wimmerer <wrkimmerer@outlook.com> and contributors
    =#
    n = 5
    A0 = sparse(SuiteSparse.increment!([0,4,1,1,2,2,0,1,2,3,4,4]),
                SuiteSparse.increment!([0,4,0,2,1,2,1,4,3,2,1,2]),
                [2.,1.,3.,4.,-1.,-3.,3.,6.,2.,1.,4.,2.], n, n)
    A1 = sparse(increment!([0,4,1,1,2,2,0,1,2,3,4,4]),
                increment!([0,4,0,2,1,2,1,4,3,2,1,2]),
                [2.,1.,3.,4.,-1.,-3.,3.,9.,2.,1.,4.,2.], n, n)
    b = [8., 45., -3., 3., 19.]
    #= End of MIT LICENSE =#

    sol1 = A0  \ b
    sol2 = A0' \ b
    sol3 = A1  \ b

    A0_d = CuSparseMatrixCSR(A0)
    A1_d = CuSparseMatrixCSR(A1)

    @testset "Normal mode ($sym)" for sym in [:KLU, :RF]
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

    @testset "Batch mode ($sym)" for sym in [:KLU, :RF]
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
        @test Matrix(B_d) ≈ X1_sol

        # Refactorization
        lu!(rf, A1_d)

        ldiv!(x_d, rf, b_d)
        @test Matrix(x_d) ≈ X2_sol
    end
end

