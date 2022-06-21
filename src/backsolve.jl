
struct CuSparseSV
    transa::CUSPARSE.SparseChar
    descL::CUSPARSE.CuMatrixDescriptor
    descU::CUSPARSE.CuMatrixDescriptor
    infoL::Vector{Ptr{Cvoid}}
    infoU::Vector{Ptr{Cvoid}}
    buffer::CuVector{UInt8}
end

function CuSparseSV(
    A::CUSPARSE.CuSparseMatrixCSR, transa::CUSPARSE.SparseChar,
)
    # Lower triangular part
    descL = CUSPARSE.CuMatrixDescriptor('G', 'L', 'U', 'O')
    # Upper triangular part
    descU = CUSPARSE.CuMatrixDescriptor('G', 'U', 'N', 'O')
    m, n = A.dims

    infoL = CUSPARSE.csrsv2Info_t[0]
    CUSPARSE.cusparseCreateCsrsv2Info(infoL)

    # Compute buffer size
    outL = Ref{Cint}(1)
    CUSPARSE.cusparseDcsrsv2_bufferSize(
        CUSPARSE.handle(), transa, m, SparseArrays.nnz(A),
        descL, SparseArrays.nonzeros(A), A.rowPtr, A.colVal, infoL[1],
        outL,
    )

    infoU = CUSPARSE.csrsv2Info_t[0]
    CUSPARSE.cusparseCreateCsrsv2Info(infoU)
    outU = Ref{Cint}(1)
    CUSPARSE.cusparseDcsrsv2_bufferSize(
        CUSPARSE.handle(), transa, m, SparseArrays.nnz(A),
        descU, SparseArrays.nonzeros(A), A.rowPtr, A.colVal, infoU[1],
        outU,
    )

    # Allocate buffer
    @assert outL[] == outU[]
    n_bytes = outL[]::Cint
    buffer = CUDA.zeros(UInt8, n_bytes)

    # Allocate triangular
    for (desc, info) in [(descL, infoL), (descU, infoU)]
        CUSPARSE.cusparseDcsrsv2_analysis(
            CUSPARSE.handle(), transa, m, SparseArrays.nnz(A),
            desc, SparseArrays.nonzeros(A), A.rowPtr, A.colVal, info[1],
            CUSPARSE.CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer,
        )
        posit = Ref{Cint}(1)
        CUSPARSE.cusparseXcsrsv2_zeroPivot(CUSPARSE.handle(), info[1], posit)
        if posit[] >= 0
            error("Structural/numerical zero in A at ($(posit[]),$(posit[])))")
        end
    end

    return CuSparseSV(transa, descL, descU, infoL, infoU, buffer)
end

function backsolve!(s::CuSparseSV, A::CUSPARSE.CuSparseMatrixCSR, X::CuVector)
    m,n = A.dims
    alpha = 1.0

    operations = if s.transa == 'N'
        [(s.descL, s.infoL), (s.descU, s.infoU)]
    elseif s.transa == 'T'
        [(s.descU, s.infoU), (s.descL, s.infoL)]
    end

    for (desc, info) in operations
        CUSPARSE.cusparseDcsrsv2_solve(
            CUSPARSE.handle(),
            s.transa, m, SparseArrays.nnz(A), alpha, desc,
            SparseArrays.nonzeros(A), A.rowPtr, A.colVal, info[1],
            X, X,
            CUSPARSE.CUSPARSE_SOLVE_POLICY_USE_LEVEL, s.buffer,
        )
    end
end


struct CuSparseSM
    transa::CUSPARSE.SparseChar
    descL::CUSPARSE.CuMatrixDescriptor
    descU::CUSPARSE.CuMatrixDescriptor
    infoL::Vector{Ptr{Cvoid}}
    infoU::Vector{Ptr{Cvoid}}
    buffer::CuVector{UInt8}
end

function CuSparseSM(
    A::CUSPARSE.CuSparseMatrixCSR, transa::CUSPARSE.SparseChar, X::CuMatrix,
)
    descL = CUSPARSE.CuMatrixDescriptor('G', 'L', 'U', 'O')
    descU = CUSPARSE.CuMatrixDescriptor('G', 'U', 'N', 'O')
    m, n = A.dims
    transxy = 'N'
    alpha = 1.0
    nX = size(X, 2)
    ldx = max(1, stride(X, 2))

    infoL = CUSPARSE.csrsm2Info_t[0]
    CUSPARSE.cusparseCreateCsrsm2Info(infoL)

    outL = Ref{UInt64}(1)
    # TODO
    CUSPARSE.cusparseDcsrsm2_bufferSizeExt(
        CUSPARSE.handle(), 0, transa, transxy, m, nX, SparseArrays.nnz(A),
        alpha, descL, SparseArrays.nonzeros(A), A.rowPtr, A.colVal, X, ldx, infoL[1],
        CUSPARSE.CUSPARSE_SOLVE_POLICY_USE_LEVEL,
        outL,
    )

    infoU = CUSPARSE.csrsm2Info_t[0]
    CUSPARSE.cusparseCreateCsrsm2Info(infoU)
    outU = Ref{UInt64}(1)
    CUSPARSE.cusparseDcsrsm2_bufferSizeExt(
        CUSPARSE.handle(), 0, transa, transxy, m, nX, SparseArrays.nnz(A),
        alpha, descU, SparseArrays.nonzeros(A), A.rowPtr, A.colVal, X, ldx, infoU[1],
        CUSPARSE.CUSPARSE_SOLVE_POLICY_USE_LEVEL,
        outU,
    )

    @assert outL[] == outU[]
    n_bytes = outL[]::UInt64
    buffer = CUDA.zeros(UInt8, n_bytes)

    for (desc, info) in [(descL, infoL), (descU, infoU)]
        CUSPARSE.cusparseDcsrsm2_analysis(
            CUSPARSE.handle(), 0, transa, transxy, m, nX, SparseArrays.nnz(A), alpha,
            desc, SparseArrays.nonzeros(A), A.rowPtr, A.colVal, X, ldx, info[1],
            CUSPARSE.CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer,
        )
        posit = Ref{Cint}(1)
        CUSPARSE.cusparseXcsrsm2_zeroPivot(CUSPARSE.handle(), info[1], posit)

        if posit[] >= 0
            error("Structural/numerical zero in A at ($(posit[]),$(posit[])))")
        end
    end

    return CuSparseSM(transa, descL, descU, infoL, infoU, buffer)
end

function backsolve!(s::CuSparseSM, A::CUSPARSE.CuSparseMatrixCSR, X::CuMatrix)
    m,n = A.dims
    alpha = 1.0
    transxy = 'N'
    nX = size(X, 2)
    ldx = max(1, stride(X, 2))

    operations = if s.transa == 'N'
        [(s.descL, s.infoL), (s.descU, s.infoU)]
    elseif s.transa == 'T'
        [(s.descU, s.infoU), (s.descL, s.infoL)]
    end

    for (desc, info) in operations
        CUSPARSE.cusparseDcsrsm2_solve(
            CUSPARSE.handle(), 0, s.transa, transxy, m, nX, SparseArrays.nnz(A), alpha,
            desc, SparseArrays.nonzeros(A), A.rowPtr, A.colVal, X, ldx, info[1],
            CUSPARSE.CUSPARSE_SOLVE_POLICY_USE_LEVEL, s.buffer,
        )
    end
end

