
struct GLULowLevel{Ti, Tv}
    glu::csrgluInfo_t
    n::Ti
    nnzA::Ti
    drowsA::CuVector{Ti}
    dcolsA::CuVector{Ti}
    dvalsA::CuVector{Tv}
    rhs::CuVector{Tv}
    dT::CuVector{UInt8}
end

function GLULowLevel(
    lu_host::RFSymbolicAnalysis{T, Ti}
) where {T, Ti}
    n, m = size(lu_host)
    glu = cusolverGluCreate()

    # Get matrix M in CSR format from symbolic factorization
    rowsM, colsM = compute_M(lu_host)
    nnzM = rowsM[n+1]

    # Load data on device
    rhs = CUDA.zeros(T, n)
    drowsA = CuVector(lu_host.rowsA)
    dcolsA = CuVector(lu_host.colsA)
    dvalsA = CuVector(lu_host.valsA)

    descA = CUSPARSE.CuMatrixDescriptor('G', 'L', 'N', 'Z')
    descM = CUSPARSE.CuMatrixDescriptor('G', 'L', 'N', 'Z')
    for desc in [descA, descM]
        CUSPARSE.cusparseSetMatType(desc, CUSPARSE.CUSPARSE_MATRIX_TYPE_GENERAL)
        CUSPARSE.cusparseSetMatIndexBase(desc, CUSPARSE.CUSPARSE_INDEX_BASE_ZERO)
    end
    # Setup GLU instance.
    cusolverSpDgluSetup(
        CUSOLVER.sparse_handle(),
        n,
        lu_host.nnzA,
        descA,
        lu_host.rowsA,
        lu_host.colsA,
        lu_host.P,
        lu_host.Q,
        nnzM,
        descM,
        rowsM,
        colsM,
        glu,
    )
    # Create buffer.
    size_buffer = Ref{Csize_t}(0)
    cusolverSpDgluBufferSize(
        CUSOLVER.sparse_handle(),
        glu,
        size_buffer,
    )
    dT = CUDA.zeros(UInt8, size_buffer[])
    # Analysis.
    cusolverSpDgluAnalysis(
        CUSOLVER.sparse_handle(),
        glu,
        dT,
    )
    # Move the factors on the GPU.
    cusolverSpDgluReset(
        CUSOLVER.sparse_handle(),
        n,
        lu_host.nnzA,
        descA,
        dvalsA,
        drowsA,
        dcolsA,
        glu,
    )
    # Factorize.
    cusolverSpDgluFactor(
        CUSOLVER.sparse_handle(),
        glu,
        dT,
    )
    return GLULowLevel(
        glu, n, lu_host.nnzA, drowsA, dcolsA, dvalsA, rhs, dT,
    )
end

function glu_refactor!(glu::GLULowLevel, A::CUSPARSE.CuSparseMatrixCSR)
    descA = CUSPARSE.CuMatrixDescriptor('G', 'L', 'N', 'Z')
    copyto!(glu.dvalsA, A.nzVal)
    cusolverSpDgluReset(
        CUSOLVER.sparse_handle(),
        glu.n,
        glu.nnzA,
        descA,
        A.nzVal,
        glu.drowsA, glu.dcolsA,
        glu.glu,
    )
    cusolverSpDgluFactor(
        CUSOLVER.sparse_handle(),
        glu.glu,
        glu.dT,
    )
    return
end

function glu_solve!(glu::GLULowLevel, x::CuVector)
    descA = CUSPARSE.CuMatrixDescriptor('G', 'L', 'N', 'Z')
    ite_refine_succ = Ref{Cint}(0)
    r_nrminf = Ref{Cdouble}(0.0)
    cusolverSpDgluSolve(
        CUSOLVER.sparse_handle(),
        glu.n,
        glu.nnzA,
        descA,
        glu.dvalsA, glu.drowsA, glu.dcolsA,
        x,  # RHS
        glu.rhs,  # LHS
        ite_refine_succ,
        r_nrminf,
        glu.glu,
        glu.dT,
    )
    copyto!(x, glu.rhs)
    return
end

