
struct RFSymbolicAnalysis{Tv, Ti}
    n::Ti
    m::Ti
    nnzA::Ti
    rowsA::Vector{Ti}
    colsA::Vector{Ti}
    valsA::Vector{Tv}
    nnzL::Ti
    rowsL::Vector{Ti}
    colsL::Vector{Ti}
    valsL::Vector{Tv}
    nnzU::Ti
    rowsU::Vector{Ti}
    colsU::Vector{Ti}
    valsU::Vector{Tv}
    P::Vector{Ti}
    Q::Vector{Ti}
end

Base.size(rf::RFSymbolicAnalysis) = (rf.n, rf.m)

# By default, we run the symbolic analysis on the CPU using
# the low-level utilities provided in cusolver.
function rf_symbolic_analysis(
    A::CUSPARSE.CuSparseMatrixCSR{T, Ti};
    ordering=:AMD, tol=1e-8, check=true,
) where {T, Ti}
    m, n = size(A)
    @assert m == n # only squared matrices are supported
    nnzA = SparseArrays.nnz(A)

    # Transfer data to host
    h_rowsA, h_colsA, h_valsA = get_csr_host(A)

    # cusolverRf is 0-based
    decrement!(h_rowsA)
    decrement!(h_colsA)
    h_Qreorder = zeros(Cint, n)
    # Create duplicate matrix for reordering
    h_rowsB = copy(h_rowsA)
    h_colsB = copy(h_colsA)
    h_valsB = copy(h_valsA)

    spH = CUSOLVER.sparse_handle()

    # Create matrix descriptor
    desca = CUSPARSE.CuMatrixDescriptor()
    CUSPARSE.cusparseSetMatType(desca, CUSPARSE.CUSPARSE_MATRIX_TYPE_GENERAL)
    CUSPARSE.cusparseSetMatIndexBase(desca, CUSPARSE.CUSPARSE_INDEX_BASE_ZERO)

    # Reordering
    if ordering == :AMD
        CUSOLVER.cusolverSpXcsrsymamdHost(
            spH,
            n, nnzA, desca,
            h_rowsA, h_colsA, h_Qreorder,
        )
    elseif ordering == :MDQ
        CUSOLVER.cusolverSpXcsrsymmdqHost(
            spH,
            n, nnzA, desca,
            h_rowsA, h_colsA, h_Qreorder,
        )
    elseif ordering == :METIS
        CUSOLVER.cusolverSpXcsrmetisndHost(
            spH,
            n, nnzA, desca,
            h_rowsA, h_colsA, C_NULL, h_Qreorder,
        )
    elseif ordering == :RCM
        CUSOLVER.cusolverSpXcsrsymrcmHost(
            spH,
            n, nnzA, desca,
            h_rowsA, h_colsA, h_Qreorder,
        )
    end

    h_mapBfromA = zeros(Cint, nnzA)
    @inbounds for i in 1:nnzA
        h_mapBfromA[i] = i # identity matrix
    end

    # Compute permutation in two steps
    size_perm = Ref{Csize_t}(0)
    CUSOLVER.cusolverSpXcsrperm_bufferSizeHost(
        spH,
        m, n, nnzA, desca,
        h_rowsB, h_colsB, h_Qreorder, h_Qreorder,
        size_perm,
    )

    buffer_cpu = zeros(Cint, size_perm[])
    CUSOLVER.cusolverSpXcsrpermHost(
        spH,
        m, n, nnzA, desca,
        h_rowsB, h_colsB, h_Qreorder, h_Qreorder, h_mapBfromA,
        buffer_cpu,
    )

    # Apply permutation
    h_valsB = h_valsA[h_mapBfromA]

    # LU Factorization
    info = Ref{CUSOLVER.csrluInfoHost_t}()
    CUSOLVER.cusolverSpCreateCsrluInfoHost(info)

    CUSOLVER.cusolverSpXcsrluAnalysisHost(
        spH,
        m, nnzA, desca,
        h_rowsB, h_colsB, info[],
    )

    size_internal = Ref{UInt64}(0)
    size_lu = Ref{UInt64}(0)
    CUSOLVER.cusolverSpDcsrluBufferInfoHost(
        spH,
        n, nnzA, desca,
        h_valsB, h_rowsB, h_colsB,
        info[],
        size_internal, size_lu,
    )

    n_bytes = size_lu[] * sizeof(Cint)
    buffer_lu = zeros(Cint, size_lu[])
    pivot_threshold = 1.0

    CUSOLVER.cusolverSpDcsrluFactorHost(
        spH, n, nnzA, desca,
        h_valsB, h_rowsB, h_colsB,
        info[], pivot_threshold,
        buffer_lu,
    )

    # Check singularity
    if check
        singularity = Ref{Cint}(0)
        CUSOLVER.cusolverSpDcsrluZeroPivotHost(
            spH, info[], tol, singularity,
        )

        # Check that the matrix is nonsingular
        if singularity[] >= 0
            throw(LinearAlgebra.SingularException(singularity[]))
        end
    end

    # Get size of L and U
    pnnzU = Ref{Cint}(0)
    pnnzL = Ref{Cint}(0)
    CUSOLVER.cusolverSpXcsrluNnzHost(
        spH,
        pnnzL, pnnzU, info[],
    )

    nnzL = pnnzL[]
    nnzU = pnnzU[]

    # Retrieve L and U matrices
    h_Plu = zeros(Cint, m)
    h_Qlu = zeros(Cint, n)

    h_valsL = zeros(nnzL)
    h_rowsL = zeros(Cint, m+1)
    h_colsL = zeros(Cint, nnzL)

    h_valsU = zeros(nnzU)
    h_rowsU = zeros(Cint, m+1)
    h_colsU = zeros(Cint, nnzU)

    # Extract
    CUSOLVER.cusolverSpDcsrluExtractHost(
        spH,
        h_Plu, h_Qlu,
        desca,
        h_valsL, h_rowsL, h_colsL,
        desca,
        h_valsU, h_rowsU, h_colsU,
        info[],
        buffer_lu,
    )

    h_P = h_Qreorder[h_Plu .+ 1]
    h_Q = h_Qreorder[h_Qlu .+ 1]

    return RFSymbolicAnalysis{T, Cint}(
        m, n, nnzA, h_rowsA, h_colsA, h_valsA,
        nnzL, h_rowsL, h_colsL, h_valsL,
        nnzU, h_rowsU, h_colsU, h_valsU,
        h_P, h_Q,
    )
end

# Evaluate L + U in 0-based indexing
function compute_M(sym::RFSymbolicAnalysis)
    n = sym.n
    nnzM = sym.nnzL + sym.nnzU

    rowsM = zeros(Cint, n + 1)
    colsM = zeros(Cint, nnzM)

    # Add L
    k = 1
    for i in 1:n
        for c in sym.rowsL[i]:sym.rowsL[i+1]-1
            j = sym.colsL[c+1]
            @assert i != j
            rowsM[i] += 1
            colsM[k] = j
            k +=1
        end
        for c in sym.rowsU[i]:sym.rowsU[i+1]-1
            j = sym.colsU[c+1]
            rowsM[i] += 1
            colsM[k] = j
            k += 1
        end
    end
    last = 1
    for i in 1:n+1
        tmp = rowsM[i]
        rowsM[i] = last
        last += tmp
    end

    # Enforce 0-based indexing
    decrement!(rowsM)

    return rowsM, colsM
end

