
using KLU

@testset "KLU symbolic factorization" begin
    n = 512
    A = sprand(n, n, .2)
    # TODO: need to add diagonal term to avoid nonzero F matrix in KLU
    A += 1e-4I #A'
    b = rand(n)
    # Compute solution with UMFPACK
    solution = A \ b

    sym_klu = CUSOLVERRF.klu_symbolic_analysis(A)
    rflu = CUSOLVERRF.RFLowLevel(sym_klu)

    d_b = CuVector{Float64}(b)
    d_x = CUDA.zeros(Float64, n)

    copyto!(d_x, d_b)
    CUSOLVERRF.rf_solve!(rflu, d_x)
    res = Array(d_x)
    @test isapprox(res, solution)

    # Test refactorization
    d_A = CuSparseMatrixCSR(A)
    scale = 2.0
    d_A.nzVal .*= scale
    CUSOLVERRF.rf_refactor!(rflu, d_A)
    copyto!(d_x, d_b)
    CUSOLVERRF.rf_solve!(rflu, d_x)
    res = Array(d_x)
    @test isapprox(res, solution ./ scale)
end

