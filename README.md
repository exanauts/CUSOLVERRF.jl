# CUSOLVERRF.jl

This package is a thin Julia wrapper for [cusolverRF](https://docs.nvidia.com/cuda/cusolver/index.html#cuSolverRF-reference`)
It extends [CUSOLVER.jl](https://github.com/JuliaGPU/CUDA.jl/tree/master/lib/cusolver), and is compatible with the rest of the Julia CUDA ecosystem.

## What is cusolverRF?

cusolverRF is a package to perform fast sparse refactorization on CUDA GPU.
The LU factorization of a sparse matrix usually happens in two steps.

1. During the symbolic factorization, the matrix is reordered to minimize the fill-in, and the sparsity patterns of the `L` and `U` factor are computed.
2. During the numerical factorization, the coefficients of the factors `L` and `U` are computed.

The first step is challenging to compute on the GPU, as most sparse matrices
have unstructured sparsity pattern. [cusolver](https://docs.nvidia.com/cuda/cusolver/index.html#cusolver-lt-t-gt-csrlsvlu) implements a default LU factorization,
but it transfers the sparse matrix on the host under the hood to perform
the symbolic factorization. This impairs the performance when the same
matrix has to be factorized multiple times.

`cusolverRF` uses a different approach: it computes the symbolic factorization on the
CPU and then transfers it on the device. If the coefficient of the sparse matrix
are updated **without affecting the sparsity pattern of the matrix**, then there
is no need to recompute the symbolic factorization and we can compute the
numerical factorization entirely on the device, **without data transfer between
the host and the device**.

Hence, `cusolverRF` is very efficient when the same sparse matrix has to be factorized
multiple times. The package is also efficient to solve a sparse linear system
with multiple right-hand-side entirely on the GPU.


## How to use cusolverRF?

Suppose we have a sparse matrix `dA` instantiated on the GPU.
```julia
using SparseArrays, CUDA, LinearAlgebra
using CUDA.CUSPARSE

# Generate a random example
n = 10
A = sprand(n, n, .2) + I 
# Move A to the GPU
dA = CuSparseMatrixCSR(A)
```
Computing the LU factorization of the sparse matrix `dA` with `cusolverRF` simply amount to
```julia
using CUSOLVERRF
rf = CUSOLVERRF.RFLU(dA; symbolic=:RF)
```
In this step, the matrix is moved back to the host to compute the
symbolic factorization. The factors `L` and `U` are then deported on
the device, and can be accessed as a matrix `M = L + U`
```julia
rf.M

```
The symbolic factorization is computed by default (`symbolic=:RF`) using `cusolver` internal
routines, which can be inefficient.
As an alternative, CUSOLVERRF.jl allows to compute the symbolic factorization
with [KLU](https://github.com/JuliaSparse/KLU.jl) (`symbolic=:KLU`). However, this second option is still experimental.

Then, computing the solution of the linear system $Ax = b$ translates to
```julia
b = rand(n)      # random RHS
db = CuVector(b) # move RHS to the device
ldiv!(rf, db)    # compute solution inplace

```

Suppose now we change the coefficients of the matrix `dA`, **without
modifying its sparsity pattern**. Then, the new system can be
solved entirely on the device as:
```julia
# update coefficients (use whichever update you want here)
dA.nzVal .*= 2.0
# update factorization on the device
lu!(rf, dA)
# resolve Ax = b
copyto!(db, b)
ldiv!(rf, b)
```

## How to solve a system with multiple right-hand-side?

Any `RFLU` instance stores different buffers to avoid unnecessary
allocations when calling `lu!` and `ldiv!`. To solve a system
with `k` multiple right-hand-side $AX=B$, one has
to instantiate a `RFLU` with the proper buffers, as:
```julia
k = 64
rf_batched = CUSOLVERRF.RFLU(dA; nrhs=k)

```
Then, solving the system $AX=B$ on the GPU simply amounts to
```julia
B = rand(n, k)
dX = CuMatrix(B)
# Compute solution inplace
ldiv!(rf_batched, dX)

```

## How to use the low-level wrapper?

Someone, one prefers to have full control on the factorization,
without any boilerplate code. CUSOLVERRF provides a direct interface
to `cusolverRF` for this purpose. In that case, a `cusolverRF` instance
is instantiated as
```julia
rf_lowlevel = CUSOLVERRF.RFLowLevel(
    dA;              # matrix to factorize
    fast_mode=true,  # fast mode is recommended
    ordering=:AMD,   # :AMD, :MDQ, :METIS or :RCM
    check=true,
    factorization_algo=CUSOLVERRF.CUSOLVERRF_FACTORIZATION_ALG0,
    triangular_algo=CUSOLVERRF.CUSOLVERRF_TRIANGULAR_SOLVE_ALG1,
)

```
Once the matrix factorized,
solving the linear system $Ax =b$ translates to
```julia
b = rand(n)
dx = CuVector(b)
CUSOLVERRF.rf_solve!(rf_lowlevel, dx)

```
And refactorizing the matrix inplace on the device:
```julia
dA.nzVal .*= 2.0
CUSOLVERRF.rf_refactor!(rf_lowlevel, dA)

```

## Funding
This package was supported by the Exascale Computing Project (17-SC-20-SC), a joint project of the U.S. Department of Energy’s Office of Science and National Nuclear Security Administration, responsible for delivering a capable exascale ecosystem, including software, applications, and hardware technology, to support the nation’s exascale computing imperative.

