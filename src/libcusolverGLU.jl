
import CUDA.CUSOLVER: libcusolver
import CUDA.CUSOLVER: cusolverStatus_t, cusolverSpHandle_t
import CUDA.CUSPARSE: cusparseMatDescr_t

mutable struct cusolverGluInfo end
const csrgluInfo_t = Ptr{cusolverGluInfo}

function cusolverSpCreateGluInfo(handle)
    @ccall libcusolver.cusolverSpCreateGluInfo(handle::Ptr{csrgluInfo_t})::cusolverStatus_t
end

function cusolverSpDestroyGluInfo(handle)
    @ccall libcusolver.cusolverSpCreateGluInfo(handle::csrgluInfo_t)::cusolverStatus_t
end

function cusolverSpDgluSetup(handle, m, nnzA, descrA, Ap, Aj, P, Q, nnzM, descrM, Mp, Mj, info)
    @ccall libcusolver.cusolverSpDgluSetup(
        handle::cusolverSpHandle_t,
        m::Cint,
        nnzA::Cint,
        descrA::cusparseMatDescr_t,
        Ap::Ptr{Cint},
        Aj::Ptr{Cint},
        P::Ptr{Cint},
        Q::Ptr{Cint},
        nnzM::Cint,
        descrM::cusparseMatDescr_t,
        Mp::Ptr{Cint},
        Mj::Ptr{Cint},
        info::csrgluInfo_t,
    )::cusolverStatus_t
end

function cusolverSpDgluBufferSize(handle, info, bufferSize)
    @ccall libcusolver.cusolverSpDgluBufferSize(
        handle::cusolverSpHandle_t,
        info::csrgluInfo_t,
        bufferSize::Ptr{Csize_t},
    )::cusolverStatus_t
end

function cusolverSpDgluAnalysis(handle, info, workspace)
    @ccall libcusolver.cusolverSpDgluAnalysis(
        handle::cusolverSpHandle_t,
        info::csrgluInfo_t,
        workspace::CuPtr{Cvoid},
    )::cusolverStatus_t
end

function cusolverSpDgluReset(handle, m, nnzA, descrA, Az, Ap, Aj, info)
    @ccall libcusolver.cusolverSpDgluReset(
        handle::cusolverSpHandle_t,
        m::Cint,
        nnzA::Cint,
        descrA::cusparseMatDescr_t,
        Az::CuPtr{Cdouble},
        Ap::CuPtr{Cint},
        Aj::CuPtr{Cint},
        info::csrgluInfo_t,
    )::cusolverStatus_t
end

function cusolverSpDgluFactor(handle, info, workspace)
    @ccall libcusolver.cusolverSpDgluFactor(
        handle::cusolverSpHandle_t,
        info::csrgluInfo_t,
        workspace::CuPtr{Cvoid},
    )::cusolverStatus_t
end

function cusolverSpDgluSolve(handle, m, nnzA, descrA, Az, Ap, Aj, b, x, ite_refine_succ, r_nrminf, info, workspace)
    @ccall libcusolver.cusolverSpDgluSolve(
        handle::cusolverSpHandle_t,
        m::Cint,
        nnzA::Cint,
        descrA::cusparseMatDescr_t,
        Az::CuPtr{Cdouble},
        Ap::CuPtr{Cint},
        Aj::CuPtr{Cint},
        b::CuPtr{Cdouble},
        x::CuPtr{Cdouble},
        ite_refine_succ::Ptr{Cint},
        r_nrminf::Ptr{Cdouble},
        info::csrgluInfo_t,
        workspace::CuPtr{Cvoid},
    )::cusolverStatus_t
end


function cusolverGluCreate()
    handle_ref = Ref{csrgluInfo_t}()
    cusolverSpCreateGluInfo(handle_ref)
    return handle_ref[]
end

function cusolverGluFree(handle)
    if handle != C_NULL
        cusolverSpDestroyGluInfo(handle)
        handle = C_NULL
    end
end

