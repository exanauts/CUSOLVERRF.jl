using Test
using CUSOLVERRF
using SuiteSparse

CUSOLVERRF.CUDA.versioninfo()
include("cusolverRF.jl")
include("klu.jl")
include("interface.jl")

