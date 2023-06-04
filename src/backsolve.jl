
function _cu_matrix_description(A::CUSPARSE.CuSparseMatrixCSR, uplo, diag, index)
    desc = CUSPARSE.CuSparseMatrixDescriptor(A, index)
    cusparse_uplo = Ref{CUSPARSE.cusparseFillMode_t}(uplo)
    cusparse_diag = Ref{CUSPARSE.cusparseDiagType_t}(diag)
    CUSPARSE.cusparseSpMatSetAttribute(desc, 'F', cusparse_uplo, Csize_t(sizeof(cusparse_uplo)))
    CUSPARSE.cusparseSpMatSetAttribute(desc, 'D', cusparse_diag, Csize_t(sizeof(cusparse_diag)))
    return desc
end

struct CuSparseSV <: AbstractBacksolve
    algo::CUSPARSE.cusparseSpSVAlg_t
    transa::CUSPARSE.SparseChar
    descL::CUSPARSE.CuSparseMatrixDescriptor
    descU::CUSPARSE.CuSparseMatrixDescriptor
    infoL::CUSPARSE.CuSparseSpSVDescriptor
    infoU::CUSPARSE.CuSparseSpSVDescriptor
    buffer::CuVector{UInt8}
end

function CuSparseSV(
    A::CUSPARSE.CuSparseMatrixCSR{T}, transa::CUSPARSE.SparseChar;
    algo=CUSPARSE.CUSPARSE_SPSV_ALG_DEFAULT,
) where T
    # CUDA 12 is required for inplace computation in cusparseSpSV
    @assert CUDA.runtime_version() >= v"12.0"
    n, m = size(A)
    @assert n == m

    # Lower triangular part
    descL = _cu_matrix_description(A, 'L', 'U', 'O')
    # Upper triangular part
    descU = _cu_matrix_description(A, 'U', 'N', 'O')

    # Dummy coefficient
    alpha = T(1.0)

    x = CUDA.zeros(T, n)
    descX = CUSPARSE.CuDenseVectorDescriptor(x)

    # Descriptor for lower-triangular SpSV operation
    spsv_L = CUSPARSE.CuSparseSpSVDescriptor()
    # Compute buffer size
    outL = Ref{Csize_t}(1)
    CUSPARSE.cusparseSpSV_bufferSize(
        CUSPARSE.handle(), transa, Ref{T}(alpha), descL, descX, descX, T, algo, spsv_L, outL,
    )

    # Descriptor for lower-triangular SpSV operation
    spsv_U = CUSPARSE.CuSparseSpSVDescriptor()
    # Compute buffer size
    outU = Ref{Csize_t}(1)
    CUSPARSE.cusparseSpSV_bufferSize(
        CUSPARSE.handle(), transa, Ref{T}(alpha), descU, descX, descX, T, algo, spsv_U, outU,
    )

    # Allocate buffer
    @assert outL[] == outU[]
    n_bytes = outL[]::Csize_t
    buffer = CUDA.zeros(UInt8, n_bytes)

    return CuSparseSV(algo, transa, descL, descU, spsv_L, spsv_U, buffer)
end

function backsolve!(s::CuSparseSV, A::CUSPARSE.CuSparseMatrixCSR{T}, X::CuVector{T}) where T
    m,n = A.dims
    alpha = T(1.0)

    descX = CUSPARSE.CuDenseVectorDescriptor(X)
    operations = if s.transa == 'N'
        [(s.descL, s.infoL), (s.descU, s.infoU)]
    elseif s.transa == 'T'
        [(s.descU, s.infoU), (s.descL, s.infoL)]
    end

    for (desc, info) in operations
        CUSPARSE.cusparseSpSV_analysis(
            CUSPARSE.handle(), s.transa, Ref{T}(alpha), desc, descX, descX, T, s.algo, info, s.buffer,
        )
        CUSPARSE.cusparseSpSV_solve(
            CUSPARSE.handle(), s.transa, Ref{T}(alpha), desc, descX, descX, T, s.algo, info,
        )
    end
end


struct CuSparseSM <: AbstractBacksolve
    algo::CUSPARSE.cusparseSpSMAlg_t
    transa::CUSPARSE.SparseChar
    descL::CUSPARSE.CuSparseMatrixDescriptor
    descU::CUSPARSE.CuSparseMatrixDescriptor
    infoL::CUSPARSE.CuSparseSpSMDescriptor
    infoU::CUSPARSE.CuSparseSpSMDescriptor
    buffer::CuVector{UInt8}
end

function CuSparseSM(
    A::CUSPARSE.CuSparseMatrixCSR{T}, transa::CUSPARSE.SparseChar, X::CuMatrix{T};
    algo=CUSPARSE.CUSPARSE_SPSM_ALG_DEFAULT,
) where T
    # CUDA 12 is required for inplace computation in cusparseSpSV
    @assert CUDA.runtime_version() >= v"12.0"

    n, m = size(A)
    @assert n == m

    # Lower triangular part
    descL = _cu_matrix_description(A, 'L', 'U', 'O')
    # Upper triangular part
    descU = _cu_matrix_description(A, 'U', 'N', 'O')

    # Dummy coefficient
    alpha = T(1.0)

    transx = 'N'
    nX = size(X, 2)
    ldx = max(1, stride(X, 2))

    descX = CUSPARSE.CuDenseMatrixDescriptor(X)

    # Descriptor for lower-triangular SpSV operation
    spsm_L = CUSPARSE.CuSparseSpSMDescriptor()
    # Compute buffer size
    outL = Ref{Csize_t}(1)
    CUSPARSE.cusparseSpSM_bufferSize(
        CUSPARSE.handle(), transa, transx, Ref{T}(alpha), descL, descX, descX, T, algo, spsm_L, outL,
    )

    # Descriptor for upper-triangular SpSV operation
    spsm_U = CUSPARSE.CuSparseSpSMDescriptor()
    # Compute buffer size
    outU = Ref{Csize_t}(1)
    CUSPARSE.cusparseSpSM_bufferSize(
        CUSPARSE.handle(), transa, transx, Ref{T}(alpha), descU, descX, descX, T, algo, spsm_U, outU,
    )

    @assert outL[] == outU[]
    n_bytes = outL[]::UInt64
    buffer = CUDA.zeros(UInt8, n_bytes)

    return CuSparseSM(algo, transa, descL, descU, spsm_L, spsm_U, buffer)
end

function backsolve!(s::CuSparseSM, A::CUSPARSE.CuSparseMatrixCSR{T}, X::CuMatrix{T}) where T
    m,n = A.dims
    alpha = T(1.0)

    descX = CUSPARSE.CuDenseMatrixDescriptor(X)
    transx = 'N'

    operations = if s.transa == 'N'
        [(s.descL, s.infoL), (s.descU, s.infoU)]
    elseif s.transa == 'T'
        [(s.descU, s.infoU), (s.descL, s.infoL)]
    end

    for (desc, info) in operations
        CUSPARSE.cusparseSpSM_analysis(
            CUSPARSE.handle(), s.transa, transx, Ref{T}(alpha), desc, descX, descX, T, s.algo, info, s.buffer,
        )
        CUSPARSE.cusparseSpSM_solve(
            CUSPARSE.handle(), s.transa, transx, Ref{T}(alpha), desc, descX, descX, T, s.algo, info,
        )
    end
end

