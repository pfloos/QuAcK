#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <stdlib.h>
#include <stdio.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include "utils.h"
#include "ph_rpa.h"

void ph_drpa_sing(int nO, int nBas, int nS, double *h_eps, double *h_ERI,
                      double *h_Omega, double *h_XpY, double *h_XmY) {

    double *d_eps = NULL;
    double *d_ERI = NULL;

    int nV = nBas - nO;

    long long nBas_long = (long long) nBas;
    long long nBas4 = nBas_long * nBas_long * nBas_long * nBas_long;

    long long nS_long = (long long) nS;
    long long nS2 = nS_long * nS_long;

    float elapsedTime;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);


    check_Cuda_Errors(cudaMalloc((void**)&d_eps, nBas * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMalloc((void**)&d_ERI, nBas4 * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);

    printf("CPU->GPU transfer..\n");
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
    check_Cuda_Errors(cudaMalloc((void**)&d_A-B, nS2 * sizeof(double)), "cudaMalloc", __FILE__, __LINE__);

    cudaEventRecord(start, 0);
    ph_dRPA_ApB_sing(nO, nV, nBas, nS, d_eps, d_ERI, d_ApB);
    ph_dRPA_AmB_sing(nO, nV, nBas, nS, d_eps, d_ERI, d_AmB);
    check_Cuda_Errors(cudaGetLastError(), "cudaGetLastError", __FILE__, __LINE__);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on A & B kernels = %f msec\n", elapsedTime);


    // free memory
    check_Cuda_Errors(cudaFree(d_eps), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_ERI), "cudaFree", __FILE__, __LINE__);


    // TODO
    // diagonalize A+B and A-B
    int *d_info = NULL;
    double *d_Omega = NULL;
    check_Cuda_Errors(cudaMalloc((void**)&d_info, sizeof(int)),
        "cudaMalloc", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMalloc((void**)&d_Omega, nS * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);

    cudaEventRecord(start, 0);
    diag_dn_dsyevd(nS, d_info, d_Omega, d_A);
    check_Cuda_Errors(cudaGetLastError(), "cudaGetLastError", __FILE__, __LINE__);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on diagonalization  = %f msec\n", elapsedTime);


    // transfer data to CPU
    cudaEventRecord(start, 0);
    //int info_gpu = 0;
    //check_Cuda_Errors(cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost),
    //    "cudaMemcpy", __FILE__, __LINE__);
    //if (info_gpu != 0) {
    //    printf("Error: diag_dn_dsyevd returned error code %d\n", info_gpu);
    //    exit(EXIT_FAILURE);
    //}
    check_Cuda_Errors(cudaMemcpy(h_XpY, d_, nS2 * sizeof(double), cudaMemcpyDeviceToHost), 
        "cudaMemcpy", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMemcpy(h_XmY, d_, nS2 * sizeof(double), cudaMemcpyDeviceToHost), 
        "cudaMemcpy", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMemcpy(h_Omega, d_Omega, nS * sizeof(double), cudaMemcpyDeviceToHost), 
        "cudaMemcpy", __FILE__, __LINE__);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Time elapsed on GPU -> CPU transfer = %f msec\n", elapsedTime);

    check_Cuda_Errors(cudaFree(d_info), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_A), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_B), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_Omega), "cudaFree", __FILE__, __LINE__);


}

