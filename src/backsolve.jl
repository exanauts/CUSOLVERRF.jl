
function _cu_matrix_description(A::CUSPARSE.CuSparseMatrixCSR, uplo, diag, index)
    desc = CUSPARSE.CuSparseMatrixDescriptor(A, index)
    cusparse_uplo = Ref{CUSPARSE.cusparseFillMode_t}(uplo)
    cusparse_diag = Ref{CUSPARSE.cusparseDiagType_t}(diag)
    CUSPARSE.cusparseSpMatSetAttribute(desc, 'F', cusparse_uplo, Csize_t(sizeof(cusparse_uplo)))
    CUSPARSE.cusparseSpMatSetAttribute(desc, 'D', cusparse_diag, Csize_t(sizeof(cusparse_diag)))
    return desc
end

# TODO: Add these constructors in CUDA.jl
mutable struct CuDenseVectorDescriptor2
    handle::CUSPARSE.cusparseDnVecDescr_t

    function CuDenseVectorDescriptor2(T::DataType, n::Int)
        desc_ref = Ref{CUSPARSE.cusparseDnVecDescr_t}()
        CUSPARSE.cusparseCreateDnVec(desc_ref, n, CU_NULL, T)
        obj = new(desc_ref[])
        finalizer(CUSPARSE.cusparseDestroyDnVec, obj)
        obj
    end
end

Base.unsafe_convert(::Type{CUSPARSE.cusparseDnVecDescr_t}, desc::CuDenseVectorDescriptor2) = desc.handle

mutable struct CuDenseMatrixDescriptor2
    handle::CUSPARSE.cusparseDnMatDescr_t

    function CuDenseMatrixDescriptor2(T::DataType, m::Int, n::Int)
        desc_ref = Ref{CUSPARSE.cusparseDnMatDescr_t}()
        CUSPARSE.cusparseCreateDnMat(desc_ref, m, n, m, CU_NULL, T, 'C')
        obj = new(desc_ref[])
        finalizer(CUSPARSE.cusparseDestroyDnMat, obj)
        obj
    end
end

Base.unsafe_convert(::Type{CUSPARSE.cusparseDnMatDescr_t}, desc::CuDenseMatrixDescriptor2) = desc.handle

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
    chktrans(transa)

    # CUDA 12 is required for inplace computation in cusparseSpSV
    @assert CUDA.runtime_version() ≥ v"12.3"
    n, m = size(A)
    @assert n == m

    # Lower triangular part
    descL = _cu_matrix_description(A, 'L', 'U', 'O')
    # Upper triangular part
    descU = _cu_matrix_description(A, 'U', 'N', 'O')

    # Dummy coefficient
    alpha = one(T)

    # Dummy descriptor
    descX = CuDenseVectorDescriptor2(T, n)

    # Descriptor for lower-triangular SpSV operation
    spsv_L = CUSPARSE.CuSparseSpSVDescriptor()
    # Compute buffer size
    outL = Ref{Csize_t}(1)
    CUSPARSE.cusparseSpSV_bufferSize(
        CUSPARSE.handle(), transa, Ref{T}(alpha), descL, descX, descX, T, algo, spsv_L, outL,
    )

    # Descriptor for upper-triangular SpSV operation
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

    # Analysis
    CUSPARSE.cusparseSpSV_analysis(
        CUSPARSE.handle(), transa, Ref{T}(alpha), descL, descX, descX, T, algo, spsv_L, buffer,
    )
    CUSPARSE.cusparseSpSV_analysis(
        CUSPARSE.handle(), transa, Ref{T}(alpha), descU, descX, descX, T, algo, spsv_U, buffer,
    )

    return CuSparseSV(algo, transa, descL, descU, spsv_L, spsv_U, buffer)
end

function backsolve!(s::CuSparseSV, A::CUSPARSE.CuSparseMatrixCSR{T}, X::CuVector{T}) where T
    m,n = A.dims
    alpha = one(T)

    descX = CUSPARSE.CuDenseVectorDescriptor(X)
    if s.transa == 'N'
        CUSPARSE.cusparseSpSV_solve(
            CUSPARSE.handle(), s.transa, Ref{T}(alpha), s.descL, descX, descX, T, s.algo, s.infoL,
        )
        CUSPARSE.cusparseSpSV_solve(
            CUSPARSE.handle(), s.transa, Ref{T}(alpha), s.descU, descX, descX, T, s.algo, s.infoU,
        )
    else
        CUSPARSE.cusparseSpSV_solve(
            CUSPARSE.handle(), s.transa, Ref{T}(alpha), s.descU, descX, descX, T, s.algo, s.infoU,
        )
        CUSPARSE.cusparseSpSV_solve(
            CUSPARSE.handle(), s.transa, Ref{T}(alpha), s.descL, descX, descX, T, s.algo, s.infoL,
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
    chktrans(transa)

    # CUDA 12 is required for inplace computation in cusparseSpSM
    @assert CUDA.runtime_version() ≥ v"12.3"

    n, m = size(A)
    @assert n == m

    # Lower triangular part
    descL = _cu_matrix_description(A, 'L', 'U', 'O')
    # Upper triangular part
    descU = _cu_matrix_description(A, 'U', 'N', 'O')

    # Dummy coefficient
    alpha = one(T)

    transx = 'N'
    nX = size(X, 2)
    ldx = max(1, stride(X, 2))

    # Dummy descriptor
    mX, nX = size(X)
    descX = CuDenseMatrixDescriptor2(T, mX, nX)

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

    CUSPARSE.cusparseSpSM_analysis(
        CUSPARSE.handle(), transa, transx, Ref{T}(alpha), descL, descX, descX, T, algo, spsm_L, buffer,
    )
    CUSPARSE.cusparseSpSM_analysis(
        CUSPARSE.handle(), transa, transx, Ref{T}(alpha), descU, descX, descX, T, algo, spsm_U, buffer,
    )

    return CuSparseSM(algo, transa, descL, descU, spsm_L, spsm_U, buffer)
end

function backsolve!(s::CuSparseSM, A::CUSPARSE.CuSparseMatrixCSR{T}, X::CuMatrix{T}) where T
    m,n = A.dims
    alpha = one(T)

    descX = CUSPARSE.CuDenseMatrixDescriptor(X)
    transx = 'N'

    if s.transa == 'N'
        CUSPARSE.cusparseSpSM_solve(
            CUSPARSE.handle(), s.transa, transx, Ref{T}(alpha), s.descL, descX, descX, T, s.algo, s.infoL,
        )
        CUSPARSE.cusparseSpSM_solve(
            CUSPARSE.handle(), s.transa, transx, Ref{T}(alpha), s.descU, descX, descX, T, s.algo, s.infoU,
        )
    else
        CUSPARSE.cusparseSpSM_solve(
            CUSPARSE.handle(), s.transa, transx, Ref{T}(alpha), s.descU, descX, descX, T, s.algo, s.infoU,
        )
        CUSPARSE.cusparseSpSM_solve(
            CUSPARSE.handle(), s.transa, transx, Ref{T}(alpha), s.descL, descX, descX, T, s.algo, s.infoL,
        )
    end
end
