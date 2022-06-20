module CUSOLVERRF

import LinearAlgebra
import SparseArrays

import CUDA
import CUDA: CuPtr, CuVector, CuMatrix, CuArray
import CUDA.CUSPARSE
import CUDA.CUSOLVER
import CUDA.CUBLAS: unsafe_batch, unsafe_strided_batch

import KLU

# Headers
# TODO: remove once the wrapper merged in CUDA.jl
include("libcusolverRf_common.jl")
include("libcusolverRf_api.jl")

include("utils.jl")
include("rf_wrapper.jl")
include("klu.jl")

include("backsolve.jl")
include("interface.jl")

end # module
