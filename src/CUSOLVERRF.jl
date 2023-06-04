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

include("utils.jl")
include("rf_wrapper.jl")
include("klu.jl")

include("backsolve.jl")
include("backsolve_deprecated.jl")
include("interface.jl")

end # module
