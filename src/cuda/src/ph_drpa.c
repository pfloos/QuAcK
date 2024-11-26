#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <stdlib.h>
#include <stdio.h>
#include <cublas_v2.h>

#include "ph_drpa.h"

int ph_drpa(int nO, int nBas, double *h_eps, double *h_ERI) {

    double *d_eps;
    double *d_ERI;

    int nBas2 = nBas * nBas;
    int nBas4 = nBas2 * nBas2;


    check_Cuda_Errors(cudaMalloc((void**)&d_eps, nO * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMalloc((void**)&d_ERI, nBas4 * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);


    check_Cuda_Errors(cudaMemcpy(d_eps, h_eps, nO * sizeof(double), cudaMemcpyHostToDevice), 
        "cudaMemcpy", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMemcpy(d_ERI, h_ERI, nBas4 * sizeof(double), cudaMemcpyHostToDevice), 
        "cudaMemcpy", __FILE__, __LINE__);

    // construct A matrix
    int nS = nO * (nBas * nO);
    double *d_A;
    check_Cuda_Errors(cudaMalloc((void**)&d_A, nS * nS * sizeof(double)), "cudaMalloc", __FILE__, __LINE__);
    phLR_dRPA_A_sing(nO, nBas, d_eps, d_ERI, d_A);
    check_Cuda_Errors(cudaGetLastError(), "cudaGetLastError", __FILE__, __LINE__);


    check_Cuda_Errors(cudaFree(d_eps), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_ERI), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_A), "cudaFree", __FILE__, __LINE__);


    return 0;
}

