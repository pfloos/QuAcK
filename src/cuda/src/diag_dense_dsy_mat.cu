#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>



extern "C" void diag_dn_dsyevd(int n, int *info, double *W, double *A) {

    cusolverDnHandle_t cusolverH = NULL;
    cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // Compute eigenvalues and eigenvectors
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER; // Upper triangular part of the matrix is stored

    int lwork = 0;
    double *work = NULL;

    //check_Cusolver_Errors(cusolverDnCreate(&cusolverH), "cusolverDnCreate", __FILE__, __LINE__);
    cusolverDnCreate(&cusolverH);

    // Query workspace size
    //check_Cusolver_Errors(cusolverDnDsyevd_bufferSize(cusolverH, jobz, uplo, n, A, n, W, &lwork),
    //    "cusolverDnDsyevd_bufferSize", __FILE__, __LINE__);
    //check_Cuda_Errors(cudaMalloc((void**)&work, sizeof(double) * lwork),
    //    "cudaMemcpy", __FILE__, __LINE__);
    cusolverDnDsyevd_bufferSize(cusolverH, jobz, uplo, n, A, n, W, &lwork);
    cudaMalloc((void**)&work, sizeof(double) * lwork);

    // Compute eigenvalues and eigenvectors
    //check_Cusolver_Errors(cusolverDnDsyevd(cusolverH, jobz, uplo, n, A, n, W, work, lwork, info),
    //    "cusolverDnDsyevd", __FILE__, __LINE__);
    cusolverDnDsyevd(cusolverH, jobz, uplo, n, A, n, W, work, lwork, info);

    // Clean up
    //check_Cuda_Errors(cudaFree(work), "cudaFree", __FILE__, __LINE__);
    //check_Cusolver_Errors(cusolverDnDestroy(cusolverH), "cusolverDnDestroy", __FILE__, __LINE__);

    cudaFree(work);
    cusolverDnDestroy(cusolverH);

}

