import CUDA.CUSPARSE: CuSparseMatrixCSR

"""
    RFLU <: LinearAlgebra.Factorization

Sparse LU factorization on the GPU, using
[cusolverRF](https://docs.nvidia.com/cuda/cusolver/index.html#cuSolverRF-reference).

"""
struct RFLU{Tv} <: LinearAlgebra.Factorization{Tv}
    rf::RFLowLevel{Tv}
    n::Int
    # M = L + U
    M::CuSparseMatrixCSR{Tv, Int32}
    # Permutation matrices
    P::CuSparseMatrixCSR{Tv, Int32}
    Q::CuSparseMatrixCSR{Tv, Int32}
    # buffers
    r::CuVector{Tv}
    T::CuMatrix{Tv}
    # Store backsolve to avoid allocating
    # a new structure for each backsolve
    tsv::CuSparseSV
    dsm::CuSparseSM
    tsm::CuSparseSM
end

Base.size(rf::RFLU) = size(rf.M)
Base.size(rf::RFLU, dim::Integer) = size(rf.M, dim)

LinearAlgebra.adjoint(rf::RFLU) = LinearAlgebra.Adjoint(rf)

function RFLU(
    A::Union{SparseArrays.SparseMatrixCSC{Tv}, CUSPARSE.CuSparseMatrixCSR{Tv}};
    nrhs=1, symbolic=:KLU, sym_options...
) where Tv <: Float64
    # Step 1: symbolic factorization
    sym_lu = if symbolic == :KLU
        klu_symbolic_analysis(A; sym_options...)
    elseif symbolic == :RF
        rf_symbolic_analysis(A; sym_options...)
    else
        error("Symbolic factorization $(symbolic) is not supported.")
    end

    # Step 2: instantiate cusolverRF
    rf = RFLowLevel(sym_lu)

    # Step 3: structure for triangular solves
    n = size(A, 1)
    # Bundled factors M = L + U
    M = rf_extract_factors(rf, n)
    # Permutation matrices
    p = Array(rf.dP) ; increment!(p)
    q = Array(rf.dQ) ; increment!(q)
    P_cpu = SparseArrays.sparse(1:n, p, ones(n), n, n)
    Q_cpu = SparseArrays.sparse(q, 1:n, ones(n), n, n)
    P = CUSPARSE.CuSparseMatrixCSR(P_cpu)
    Q = CUSPARSE.CuSparseMatrixCSR(Q_cpu)

    # Buffers
    r = CuVector{Tv}(undef, n)       ; fill!(r, zero(Tv))
    T = CuMatrix{Tv}(undef, n, nrhs) ; fill!(T, zero(Tv))

    tsv = CuSparseSV(M, 'T')
    dsm = CuSparseSM(M, 'N', T)
    tsm = CuSparseSM(M, 'T', T)

    return RFLU(rf, n, M, P, Q, r, T, tsv, dsm, tsm)
end

# Refactoring
function LinearAlgebra.lu!(rf::RFLU, J::CuSparseMatrixCSR)
    rf_refactor!(rf.rf, J)
end

# Direct solve
function LinearAlgebra.ldiv!(
    rf::RFLU, x::AbstractVector,
)
    rf_solve!(rf.rf, x)
end

function LinearAlgebra.ldiv!(
    y::AbstractVector, rf::RFLU, x::AbstractVector,
)
    copyto!(y, x)
    rf_solve!(rf.rf, y)
end

function LinearAlgebra.ldiv!(
    Y::AbstractMatrix, rf::RFLU, X::AbstractMatrix,
)
    @assert size(Y, 2) == size(X, 2) == size(rf.T, 2)
    Z = rf.T
    LinearAlgebra.mul!(Z, rf.P, X)
    backsolve!(rf.dsm, rf.M, Z)
    LinearAlgebra.mul!(Y, rf.Q, Z)
end

function LinearAlgebra.ldiv!(
    rf::RFLU, X::AbstractMatrix,
)
    @assert size(X, 2) == size(rf.T, 2)
    Z = rf.T
    LinearAlgebra.mul!(Z, rf.P, X)
    backsolve!(rf.dsm, rf.M, Z)
    LinearAlgebra.mul!(X, rf.Q, Z)
end

# Backward solve
function LinearAlgebra.ldiv!(
    arf::LinearAlgebra.Adjoint{T, RFLU{T}}, x::AbstractVector{T},
) where T
    rf = arf.parent
    z = rf.r
    LinearAlgebra.mul!(z, rf.Q', x)
    backsolve!(rf.tsv, rf.M, z)
    LinearAlgebra.mul!(x, rf.P', z)
end

function LinearAlgebra.ldiv!(
    y::AbstractVector{T}, arf::LinearAlgebra.Adjoint{T, RFLU{T}}, x::AbstractVector{T},
) where T
    rf = arf.parent
    z = rf.r
    LinearAlgebra.mul!(z, rf.Q', x)
    backsolve!(rf.tsv, rf.M, z)
    LinearAlgebra.mul!(y, rf.P', z)
end

function LinearAlgebra.ldiv!(
    Y::AbstractMatrix{T}, arf::LinearAlgebra.Adjoint{T, RFLU{T}}, X::AbstractMatrix{T},
) where T
    rf = arf.parent
    @assert size(Y, 2) == size(X, 2) == size(rf.T, 2)
    Z = rf.T
    LinearAlgebra.mul!(Z, rf.Q', X)
    backsolve!(rf.tsm, rf.M, Z)
    LinearAlgebra.mul!(Y, rf.P', Z)
end

function LinearAlgebra.ldiv!(
    arf::LinearAlgebra.Adjoint{T, RFLU{T}}, X::AbstractMatrix{T},
) where T
    rf = arf.parent
    @assert size(X, 2) == size(rf.T, 2)
    Z = rf.T
    LinearAlgebra.mul!(Z, rf.Q', X)
    backsolve!(rf.tsm, rf.M, Z)
    LinearAlgebra.mul!(X, rf.P', Z)
end

