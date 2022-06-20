using LinearAlgebra, SparseArrays, Test
using CUDA
using CUDA.CUSOLVER
using CUDA.CUSPARSE
using CUSOLVERRF

n = 512

@testset "cusolverRF" begin
    @testset "RFLowLevel factorization" begin
        A = sprand(n, n, .2)
        A += A'
        b = rand(n)
        # Compute solution with UMFPACK
        solution = A \ b

        d_A = CuSparseMatrixCSR(A)
        d_b = CuVector{Float64}(b)
        d_x = CUDA.zeros(Float64, n)

        rflu = CUSOLVERRF.RFLowLevel(d_A)

        copyto!(d_x, d_b)
        CUSOLVERRF.rf_solve!(rflu, d_x)
        res = Array(d_x)
        @test isapprox(res, solution)
        # Test refactoring
        scale = 2.0
        d_A.nzVal .*= scale
        CUSOLVERRF.rf_refactor!(rflu, d_A)
        copyto!(d_x, d_b)
        CUSOLVERRF.rf_solve!(rflu, d_x)
        res = Array(d_x)
        @test isapprox(res, solution ./ scale)
    end

    @testset "RFBacthedLowLevel factorization" begin
        # One matrix, multiple RHS
        A = sprand(n, n, .2)
        A += A'
        nbatch = 32
        B = rand(n, nbatch)
        # Compute solution with UMFPACK
        solution = A \ B

        d_A = CuSparseMatrixCSR(A)
        d_B = CuMatrix{Float64}(B)
        d_X = CUDA.zeros(Float64, n, nbatch)
        rflu = CUSOLVERRF.RFBacthedLowLevel(d_A, nbatch)

        copyto!(d_X, d_B)
        CUSOLVER.rf_batch_solve!(rflu, d_X)
        res = Array(d_X)
        @test isapprox(res, solution)

        # Refactoring
        scale = 2.0
        d_A.nzVal .*= scale
        CUSOLVERRF.rf_batch_refactor!(rflu, d_A)
        copyto!(d_X, d_B)
        CUSOLVERRF.rf_batch_solve!(rflu, d_X)
        res = Array(d_X)
        @test isapprox(res, solution ./ scale)

        # Matrices should have the same sparsity pattern
        I, J, V = findnz(A)
        nnzA = length(V)
        # Create a batch of matrices
        As_batch = [sparse(I, J, randn(nnzA)) for i in 1:nbatch]
        d_As_batch = [CuSparseMatrixCSR(Ab) for Ab in As_batch]
        CUSOLVERRF.rf_batch_refactor!(rflu, d_As_batch)
        # ldiv!(d_X, rflu, d_B)
        # res = Array(d_X)
        # for i in 1:nbatch
        #     solution = As_batch[i] \ B[:, i]
        #     @test isapprox(res[:, i], solution)
        # end
    end
end

