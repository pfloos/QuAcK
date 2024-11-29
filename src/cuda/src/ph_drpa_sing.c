#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <stdlib.h>
#include <stdio.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include "utils.h"
#include "ph_rpa.h"
#include "my_linalg.h"





void ph_drpa_sing(int nO, int nBas, int nS, double *h_eps, double *h_ERI,
                  double *h_Omega, double *h_XpY, double *h_XmY) {


    double *d_eps = NULL;
    double *d_ERI = NULL;

    int nV = nBas - nO;

    long long nBas_long = (long long) nBas;
    long long nBas4 = nBas_long * nBas_long * nBas_long * nBas_long;

    long long nS_long = (long long) nS;
    long long nS2 = nS_long * nS_long;


    cublasHandle_t handle;
    const double alpha=1.0, beta=0.0;


    float elapsedTime;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);



    check_Cuda_Errors(cudaMalloc((void**)&d_eps, nBas * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMalloc((void**)&d_ERI, nBas4 * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);

    cudaEventRecord(start, 0);
    check_Cuda_Errors(cudaMemcpy(d_eps, h_eps, nBas * sizeof(double), cudaMemcpyHostToDevice), 
        "cudaMemcpy", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMemcpy(d_ERI, h_ERI, nBas4 * sizeof(double), cudaMemcpyHostToDevice), 
        "cudaMemcpy", __FILE__, __LINE__);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on CPU->GPU transfer = %f msec\n", elapsedTime);

    // construct A+B & A-B
    double *d_ApB = NULL;
    double *d_AmB = NULL;
    check_Cuda_Errors(cudaMalloc((void**)&d_ApB, nS2 * sizeof(double)), "cudaMalloc", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMalloc((void**)&d_AmB, nS2 * sizeof(double)), "cudaMalloc", __FILE__, __LINE__);

    cudaEventRecord(start, 0);
    ph_dRPA_ApB_sing(nO, nV, nBas, nS, d_eps, d_ERI, d_ApB);
    ph_dRPA_AmB_sing(nO, nV, nBas, nS, d_eps, d_ERI, d_AmB);
    check_Cuda_Errors(cudaGetLastError(), "cudaGetLastError", __FILE__, __LINE__);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on AmB & ApB = %f msec\n", elapsedTime);


    // free memory
    check_Cuda_Errors(cudaFree(d_eps), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_ERI), "cudaFree", __FILE__, __LINE__);


    // diagonalize A-B
    int *d_info1 = NULL;
    check_Cuda_Errors(cudaMalloc((void**)&d_info1, sizeof(int)), "cudaMalloc", __FILE__, __LINE__);
    double *d_Omega = NULL;
    check_Cuda_Errors(cudaMalloc((void**)&d_Omega, nS * sizeof(double)), "cudaMalloc", __FILE__, __LINE__);
    cudaEventRecord(start, 0);
    diag_dn_dsyevd(nS, d_info1, d_Omega, d_AmB);
    check_Cuda_Errors(cudaGetLastError(), "cudaGetLastError", __FILE__, __LINE__);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on diag AmB = %f msec\n", elapsedTime);


    // d_Omega <-- d_Omega^{0.5}
    // TODO: nb of <= 0 elements
    cudaEventRecord(start, 0);
    elementwise_dsqrt_inplace(nS, d_Omega);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on elementwise_dsqrt_inplace %f msec\n", elapsedTime);


    // d_AmBSq    = d_AmB (d_Omega)^{+0.5} (d_AmB)^T
    // d_AmBSqInv = d_AmB (d_Omega)^{-0.5} (d_AmB)^T
    cudaEventRecord(start, 0);
    double *d_AmBSq = NULL;
    check_Cuda_Errors(cudaMalloc((void**)&d_AmBSq, nS * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);
    double *d_AmBSqInv = NULL;
    check_Cuda_Errors(cudaMalloc((void**)&d_AmBSqInv, nS * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);
    A_D_At(nS, d_AmB, d_Omega, d_AmBSq);
    A_Dinv_At(nS, d_AmB, d_Omega, d_AmBSqInv);
    check_Cuda_Errors(cudaGetLastError(), "cudaGetLastError", __FILE__, __LINE__);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on d_AmBSq & d_AmBSqInv = %f msec\n", elapsedTime);




    // Dgemm
    cudaEventRecord(start, 0);
    check_Cublas_Errors(cublasCreate(&handle), "cublasCreate", __FILE__, __LINE__);

    // X + Y
    check_Cublas_Errors(cublasDgemm(handle,
                                    CUBLAS_OP_N, CUBLAS_OP_N,
                                    nS, nS, nS,
                                    &alpha,
                                    d_ApB, nS,
                                    d_AmBSq, nS,
                                    &beta,
                                    d_AmB, nS),
        "cublasDgemm", __FILE__, __LINE__);

    check_Cuda_Errors(cudaDeviceSynchronize(), "cudaDeviceSynchronize", __FILE__, __LINE__);

    // X - Y
    check_Cublas_Errors(cublasDgemm(handle,
                                    CUBLAS_OP_N, CUBLAS_OP_N,
                                    nS, nS, nS,
                                    &alpha,
                                    d_AmBSq, nS,
                                    d_AmB, nS,
                                    &beta,
                                    d_ApB, nS),
        "cublasDgemm", __FILE__, __LINE__);

    check_Cublas_Errors(cublasDestroy(handle), "cublasDestroy", __FILE__, __LINE__);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on cublasDgemm = %f msec\n", elapsedTime);




    // diagonalize
    int *d_info2 = NULL;
    check_Cuda_Errors(cudaMalloc((void**)&d_info2, sizeof(int)), "cudaMalloc", __FILE__, __LINE__);
    cudaEventRecord(start, 0);
    diag_dn_dsyevd(nS, d_info2, d_Omega, d_ApB);
    check_Cuda_Errors(cudaGetLastError(), "cudaGetLastError", __FILE__, __LINE__);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on diag ApB = %f msec\n", elapsedTime);




    // d_Omega <-- d_Omega^{0.5}
    // TODO: nb of <= 0 elements
    cudaEventRecord(start, 0);
    elementwise_dsqrt_inplace(nS, d_Omega);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on elementwise_dsqrt_inplace %f msec\n", elapsedTime);





    // Dgemm
    cudaEventRecord(start, 0);
    check_Cublas_Errors(cublasCreate(&handle), "cublasCreate", __FILE__, __LINE__);

    // X + Y
    check_Cublas_Errors(cublasDgemm(handle,
                                    CUBLAS_OP_T, CUBLAS_OP_N,
                                    nS, nS, nS,
                                    &alpha,
                                    d_ApB, nS,
                                    d_AmBSq, nS,
                                    &beta,
                                    d_AmB, nS),
        "cublasDgemm", __FILE__, __LINE__);

    check_Cuda_Errors(cudaDeviceSynchronize(), "cudaDeviceSynchronize", __FILE__, __LINE__);

    // X - Y
    check_Cublas_Errors(cublasDgemm(handle,
                                    CUBLAS_OP_T, CUBLAS_OP_N,
                                    nS, nS, nS,
                                    &alpha,
                                    d_ApB, nS,
                                    d_AmBSqInv, nS,
                                    &beta,
                                    d_AmBSq, nS),
        "cublasDgemm", __FILE__, __LINE__);

    check_Cublas_Errors(cublasDestroy(handle), "cublasDestroy", __FILE__, __LINE__);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on cublasDgemm = %f msec\n", elapsedTime);



    cudaEventRecord(start, 0);
    elementwise_dsqrt(nS, d_Omega, d_AmBSq); // avoid addition memory allocation
    A_Dinv_inplace(nS, d_AmB, d_AmBSq); // X + Y
    A_D_inplace(nS, d_ApB, d_AmBSq); // X - Y
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on final X+Y and X-Y trans = %f msec\n", elapsedTime);





    // transfer data to CPU
    cudaEventRecord(start, 0);
    check_Cuda_Errors(cudaMemcpy(h_XpY, d_AmB, nS2 * sizeof(double), cudaMemcpyDeviceToHost), 
        "cudaMemcpy", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMemcpy(h_XmY, d_ApB, nS2 * sizeof(double), cudaMemcpyDeviceToHost), 
        "cudaMemcpy", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMemcpy(h_Omega, d_Omega, nS * sizeof(double), cudaMemcpyDeviceToHost), 
        "cudaMemcpy", __FILE__, __LINE__);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on GPU -> CPU transfer = %f msec\n", elapsedTime);


    check_Cuda_Errors(cudaFree(d_info1), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_info2), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_ApB), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_AmB), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_AmBSq), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_AmBSqInv), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_Omega), "cudaFree", __FILE__, __LINE__);


}

