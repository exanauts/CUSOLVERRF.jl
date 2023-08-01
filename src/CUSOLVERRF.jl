module CUSOLVERRF

import LinearAlgebra
import SparseArrays

import CUDA
import CUDA: CuPtr, CuVector, CuMatrix, CuArray
import CUDA.CUSPARSE
import CUDA.CUSOLVER
import CUDA.CUBLAS: unsafe_batch, unsafe_strided_batch

import KLU

abstract type AbstractBacksolve end

# Core library
include("libcusolverGLU.jl")

# Low-level wrappers
include("utils.jl")
include("rf_symbolic.jl")
include("klu_symbolic.jl")
include("rf_wrapper.jl")
include("glu_wrapper.jl")

include("backsolve.jl")
include("backsolve_deprecated.jl")
include("interface.jl")

end # module
