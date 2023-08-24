
@testset "cusolverRF" begin
    # Generate data
    n = 512
    A = sprand(n, n, .2)
    b = rand(n)
    # Compute solution with UMFPACK
    solution = A \ b
    @testset "RFLowLevel factorization" begin
        d_b = CuVector{Float64}(b)
        d_x = CUDA.zeros(Float64, n)
        d_A = CuSparseMatrixCSR(A)

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
        d_A = CuSparseMatrixCSR(A)
        # One matrix, multiple RHS
        nbatch = 32
        B = rand(n, nbatch)
        # Compute solution with UMFPACK
        solution = A \ B

        d_B = CuMatrix{Float64}(B)
        d_X = CUDA.zeros(Float64, n, nbatch)
        rflu = CUSOLVERRF.RFBatchedLowLevel(d_A, nbatch)

        copyto!(d_X, d_B)
        CUSOLVERRF.rf_batch_solve!(rflu, d_X)
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
        copyto!(d_X, d_B)
        CUSOLVERRF.rf_batch_solve!(rflu, d_X)
        res = Array(d_X)
        solution = similar(res)
        for i in 1:nbatch
            solution[:, i] = As_batch[i] \ B[:, i]
        end
        @test isapprox(res, solution)
    end
end

