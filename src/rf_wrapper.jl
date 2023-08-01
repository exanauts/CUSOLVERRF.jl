
function cusolverRfCreate()
    handle_ref = Ref{CUSOLVER.cusolverRfHandle_t}()
    CUSOLVER.cusolverRfCreate(handle_ref)
    return handle_ref[]
end

function cusolverRfFree(handle)
    if handle != C_NULL
        CUSOLVER.cusolverRfDestroy(handle)
        handle = C_NULL
    end
end

mutable struct RfHandle
    handle::CUSOLVER.cusolverRfHandle_t
end

function sparse_rf_handle(;
    fast_mode=true, nzero=0.0, nboost=0.0,
    factorization_algo=CUSOLVER.CUSOLVERRF_FACTORIZATION_ALG0,
    triangular_algo=CUSOLVER.CUSOLVERRF_TRIANGULAR_SOLVE_ALG1,
)
    # Create handle
    gH = cusolverRfCreate()
    if fast_mode
        CUSOLVER.cusolverRfSetResetValuesFastMode(gH, CUSOLVER.CUSOLVERRF_RESET_VALUES_FAST_MODE_ON)
    else
        CUSOLVER.cusolverRfSetResetValuesFastMode(gH, CUSOLVER.CUSOLVERRF_RESET_VALUES_FAST_MODE_OFF)
    end
    CUSOLVER.cusolverRfSetNumericProperties(gH, nzero, nboost)
    CUSOLVER.cusolverRfSetMatrixFormat(
        gH,
        CUSOLVER.CUSOLVERRF_MATRIX_FORMAT_CSR,
        CUSOLVER.CUSOLVERRF_UNIT_DIAGONAL_ASSUMED_L
    )
    CUSOLVER.cusolverRfSetAlgs(
        gH, factorization_algo, triangular_algo,
    )
    handle = RfHandle(gH)
    finalizer(rf_free!, handle)
    return handle
end

rf_free!(rf::RfHandle) = cusolverRfFree(rf.handle)

Base.unsafe_convert(::Type{CUSOLVER.cusolverRfHandle_t}, rf::RfHandle) = rf.handle

struct RFLowLevel{T}
    rf::RfHandle
    nrhs::Int
    n::Int
    m::Int
    nnzA::Int
    drowsA::CuVector{Cint}
    dcolsA::CuVector{Cint}
    dP::CuVector{Cint}
    dQ::CuVector{Cint}
    dT::CuVector{T}
end

Base.size(rf::RFLowLevel) = (rf.n, rf.m)

# Default fallback: compute symbolic factorization with cusolver
function RFLowLevel(
    A::CUSPARSE.CuSparseMatrixCSR{T, Ti};
    ordering=:AMD, check=true, options...
) where {T, Ti}
    lu_host = rf_symbolic_analysis(A; ordering=ordering, check=check)
    return RFLowLevel(lu_host; options...)
end

function RFLowLevel(
    lu_host::RFSymbolicAnalysis{T, Ti};
    fast_mode=true,
    factorization_algo=CUSOLVER.CUSOLVERRF_FACTORIZATION_ALG0,
    triangular_algo=CUSOLVER.CUSOLVERRF_TRIANGULAR_SOLVE_ALG1,
) where {T, Ti}
    # Currently CusolverRF supports only one right-hand side.
    nrhs = 1
    n, m = size(lu_host)

    # Allocations (device)
    d_T = CUDA.zeros(Cdouble, m * nrhs)

    rf = sparse_rf_handle(;
        fast_mode=fast_mode,
        factorization_algo=factorization_algo,
        triangular_algo=triangular_algo,
    )

    # Assemble internal data structures
    CUSOLVER.cusolverRfSetupHost(
        n, lu_host.nnzA, lu_host.rowsA, lu_host.colsA, lu_host.valsA,
        lu_host.nnzL, lu_host.rowsL, lu_host.colsL, lu_host.valsL,
        lu_host.nnzU, lu_host.rowsU, lu_host.colsU, lu_host.valsU,
        lu_host.P, lu_host.Q,
        rf
    )
    # Analyze available parallelism
    CUSOLVER.cusolverRfAnalyze(rf)
    # LU factorization
    CUSOLVER.cusolverRfRefactor(rf)

    return RFLowLevel{T}(
        rf, nrhs, n, m, lu_host.nnzA,
        lu_host.rowsA, lu_host.colsA, lu_host.P, lu_host.Q, d_T
    )
end

# Update factorization inplace
function rf_refactor!(rflu::RFLowLevel{T}, A::CUSPARSE.CuSparseMatrixCSR{T, Ti}) where {T, Ti}
    CUSOLVER.cusolverRfResetValues(
        rflu.n, rflu.nnzA,
        rflu.drowsA, rflu.dcolsA, A.nzVal, rflu.dP, rflu.dQ,
        rflu.rf
    )
    CUSOLVER.cusolverRfRefactor(rflu.rf)
    return
end

# Solve system Ax = b
function rf_solve!(rflu, x::CuVector)
    n = rflu.n
    CUSOLVER.cusolverRfSolve(rflu.rf, rflu.dP, rflu.dQ, rflu.nrhs, rflu.dT, n, x, n)
    return
end

function rf_extract_factors_host(rflu::RFLowLevel, n)
    pMp = Ptr{Cint}[Ptr{Cint}(0)]
    pMj = Ptr{Cint}[Ptr{Cint}(0)]
    pMx = Ptr{Float64}[Ptr{Float64}(0)]
    pnnzM = Ref{Cint}(0)
    CUSOLVER.cusolverRfExtractBundledFactorsHost(
        rflu.rf, pnnzMM, pMp, pMj, pMx
    )
    nnzM = pnnzMM[]
    Mp = unsafe_wrap(Vector{Cint}, pMp[1], n+1)
    Mj = unsafe_wrap(Vector{Cint}, pMj[1], nnzM)
    Mx = unsafe_wrap(Vector{Float64}, pMx[1], nnzM)
    # Julia is 1-indexed
    Mp .+= Cint(1)
    Mj .+= Cint(1)
    return SparseMatrixCSC(n, n, Mp, Mj, Mx)
end

#=
    N.B.: the function wrapped in CUSOLVER.jl does not have the
    correct signature (we should use Ptr{CuPtr} instead of CuPtr{Ptr}).
    We use this custom wrapper before a fix has been ported in CUSOLVER.jl.
=#
function _cusolverRfAccessBundledFactorsDevice(handle, nnzM, Mp, Mi, Mx)
    CUSOLVER.initialize_context()
    @ccall CUSOLVER.libcusolver.cusolverRfAccessBundledFactorsDevice(handle::CUSOLVER.cusolverRfHandle_t,
                                                            nnzM::Ptr{Cint},
                                                            Mp::Ptr{CuPtr{Cint}},
                                                            Mi::Ptr{CuPtr{Cint}},
                                                            Mx::Ptr{CuPtr{Cdouble}})::CUSOLVER.cusolverStatus_t
end

function rf_extract_factors(rflu::RFLowLevel, n)
    pMp = CuPtr{Cint}[CuPtr{Cint}(0)]
    pMj = CuPtr{Cint}[CuPtr{Cint}(0)]
    pMx = CuPtr{Float64}[CuPtr{Float64}(0)]
    pnnzM = Ref{Cint}(0)
    _cusolverRfAccessBundledFactorsDevice(
        rflu.rf, pnnzM, pMp, pMj, pMx
    )
    nnzM = Int(pnnzM[])
    Mp = unsafe_wrap(CuVector{Cint}, pMp[1], n+1)
    Mj = unsafe_wrap(CuVector{Cint}, pMj[1], nnzM)
    Mx = unsafe_wrap(CuVector{Float64}, pMx[1], nnzM)
    # Avoid side effect by copying the indexings
    myMp = copy(Mp)
    myMj = copy(Mj)
    # Julia is 1-indexed
    myMp .+= Cint(1)
    myMj .+= Cint(1)
    return CuSparseMatrixCSR(myMp, myMj, Mx, (n, n))
end

# Batch factorization should not mix with classical LU factorization.
# We implement a structure apart.
struct RFBatchedLowLevel{T}
    rf::RfHandle
    batchsize::Int
    n::Int
    m::Int
    nnzA::Int
    drowsA::CuVector{Cint}
    dcolsA::CuVector{Cint}
    dP::CuVector{Cint}
    dQ::CuVector{Cint}
    dT::CuVector{T}
end

Base.size(rf::RFBatchedLowLevel) = (rf.n, rf.m)

function RFBatchedLowLevel(
    A::CUSPARSE.CuSparseMatrixCSR{T, Ti}, batchsize;
    ordering=:AMD, check=true, options...
) where {T, Ti}
    lu_host = rf_symbolic_analysis(A; ordering=ordering, check=check)
    return RFBatchedLowLevel(lu_host, batchsize; options...)
end

function RFBatchedLowLevel(
    lu_host::RFSymbolicAnalysis{T, Ti}, batchsize::Int;
    fast_mode=true,
    factorization_algo=CUSOLVER.CUSOLVERRF_FACTORIZATION_ALG0,
    triangular_algo=CUSOLVER.CUSOLVERRF_TRIANGULAR_SOLVE_ALG1,
) where {T, Ti}
    n, m = size(lu_host)

    # Allocations (device)
    d_T = CUDA.zeros(Cdouble, m * batchsize * 2)

    rf = sparse_rf_handle(;
        fast_mode=fast_mode,
        factorization_algo=factorization_algo,
        triangular_algo=triangular_algo,
    )

    # Assemble internal data structures
    h_valsA_batch = Vector{Float64}[lu_host.valsA for i in 1:batchsize]
    ptrA_batch = pointer.(h_valsA_batch)
    CUSOLVER.cusolverRfBatchSetupHost(
        batchsize,
        n, lu_host.nnzA, lu_host.rowsA, lu_host.colsA, ptrA_batch,
        lu_host.nnzL, lu_host.rowsL, lu_host.colsL, lu_host.valsL,
        lu_host.nnzU, lu_host.rowsU, lu_host.colsU, lu_host.valsU,
        lu_host.P, lu_host.Q,
        rf,
    )
    # Analyze available parallelism
    CUSOLVER.cusolverRfBatchAnalyze(rf)
    # LU factorization
    CUSOLVER.cusolverRfBatchRefactor(rf)

    return RFBatchedLowLevel{T}(
        rf, batchsize, n, m, lu_host.nnzA,
        lu_host.rowsA, lu_host.colsA, lu_host.P, lu_host.Q, d_T
    )
end

# Update factorization inplace
## Single matrix
function rf_batch_refactor!(rflu::RFBatchedLowLevel{T}, A::CUSPARSE.CuSparseMatrixCSR{T, Ti}) where {T, Ti}
    ptrs = [pointer(A.nzVal) for i in 1:rflu.batchsize]
    Aptrs = CuArray(ptrs)
    CUSOLVER.cusolverRfBatchResetValues(
        rflu.batchsize, rflu.n, rflu.nnzA,
        rflu.drowsA, rflu.dcolsA, Aptrs, rflu.dP, rflu.dQ,
        rflu.rf
    )
    CUDA.unsafe_free!(Aptrs)
    CUSOLVER.cusolverRfBatchRefactor(rflu.rf)
    return
end

## Multiple matrices
function rf_batch_refactor!(rflu::RFBatchedLowLevel{T}, As::Vector{CUSPARSE.CuSparseMatrixCSR{T, Ti}}) where {T, Ti}
    @assert length(As) == rflu.batchsize
    ptrs = [pointer(A.nzVal) for A in As]
    Aptrs = CuArray(ptrs)
    CUSOLVER.cusolverRfBatchResetValues(
        rflu.batchsize, rflu.n, rflu.nnzA,
        rflu.drowsA, rflu.dcolsA, Aptrs, rflu.dP, rflu.dQ,
        rflu.rf
    )
    CUDA.unsafe_free!(Aptrs)
    CUSOLVER.cusolverRfBatchRefactor(rflu.rf)
    return
end

function rf_batch_solve!(rflu::RFBatchedLowLevel{T}, xs::Vector{CuVector{T}}) where T
    @assert length(xs) == rflu.batchsize
    n, nrhs = rflu.n, 1
    Xptrs = unsafe_batch(xs)
    CUSOLVER.cusolverRfBatchSolve(rflu.rf, rflu.dP, rflu.dQ, nrhs, rflu.dT, n, Xptrs, n)
    CUDA.unsafe_free!(Xptrs)
    return
end

function rf_batch_solve!(rflu::RFBatchedLowLevel{T}, X::CuMatrix{T}) where T
    @assert size(X, 2) == rflu.batchsize
    n = rflu.n
    nrhs = 1
    Xptrs = unsafe_strided_batch(X)
    # Forward and backward solve
    CUSOLVER.cusolverRfBatchSolve(rflu.rf, rflu.dP, rflu.dQ, nrhs, rflu.dT, n, Xptrs, n)
    CUDA.unsafe_free!(Xptrs)
    return
end

