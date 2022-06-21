#=
    A KLU wrapper to compute the symbolic
    factorization with KLU.

=#

function klu_symbolic_analysis(
    A::SparseArrays.SparseMatrixCSC{Tv, Ti},
) where {Tv, Ti}
    n, m = size(A)

    # Initial factorization
    K = KLU.klu(A)
    K.common.scale = 0
    K.common.btf = 0
    K.common.ordering = 0
    K.common.tol = 1e-2
    # Recompute factorization with proper options
    KLU.klu!(K, A)

    if !iszero(K.F)
        error("Nonzero off-diagonal block `F` detected in factorization")
    end

    nnzA = SparseArrays.nnz(A)
    rowsA, colsA, valsA = convert2csr(A)
    rowsL, colsL, valsL = convert2csr(K.L)
    rowsU, colsU, valsU = convert2csr(K.U)

    rowsL, colsL, valsL = drop_diag_csr(rowsL, colsL, valsL)
    nnzL = length(colsL)
    nnzU = length(colsU)
    # Convert to Cint
    P = convert(Vector{Cint}, K.p)
    Q = convert(Vector{Cint}, K.q)
    rowsA = convert(Vector{Cint}, rowsA)
    colsA = convert(Vector{Cint}, colsA)
    rowsL = convert(Vector{Cint}, rowsL)
    colsL = convert(Vector{Cint}, colsL)
    rowsU = convert(Vector{Cint}, rowsU)
    colsU = convert(Vector{Cint}, colsU)

    # CUSOLVERRF is 0-based
    for vals in [rowsA, colsA, rowsL, colsL, rowsU, colsU, P, Q]
        decrement!(vals)
    end

    return RFSymbolicAnalysis{Tv, Cint}(
        n, m, nnzA, rowsA, colsA, valsA,
        nnzL, rowsL, colsL, valsL,
        nnzU, rowsU, colsU, valsU,
        P, Q,
    )
end

# Forgiving function
klu_symbolic_analysis(A::CUSPARSE.CuSparseMatrixCSR) = klu_symbolic_analysis(SparseArrays.SparseMatrixCSC(A))

