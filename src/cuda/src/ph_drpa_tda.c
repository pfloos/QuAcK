#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <stdlib.h>
#include <stdio.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include "utils.h"
#include "ph_rpa.h"

void ph_drpa_tda(int nO, int nBas, int nS, double *h_eps, double *h_ERI,
                 double *h_Omega, double *h_X) {

    double *d_eps = NULL;
    double *d_ERI = NULL;

    int nV = nBas - nO;

    int nBas2 = nBas * nBas;
    int nBas4 = nBas2 * nBas2;


    check_Cuda_Errors(cudaMalloc((void**)&d_eps, nBas * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMalloc((void**)&d_ERI, nBas4 * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);

    check_Cuda_Errors(cudaMemcpy(d_eps, h_eps, nBas * sizeof(double), cudaMemcpyHostToDevice), 
        "cudaMemcpy", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMemcpy(d_ERI, h_ERI, nBas4 * sizeof(double), cudaMemcpyHostToDevice), 
        "cudaMemcpy", __FILE__, __LINE__);

    // construct A
    double *d_A = NULL;
    check_Cuda_Errors(cudaMalloc((void**)&d_A, nS * nS * sizeof(double)), "cudaMalloc", __FILE__, __LINE__);

    ph_dRPA_A_sing(nO, nV, nBas, nS, d_eps, d_ERI, d_A);
    check_Cuda_Errors(cudaGetLastError(), "cudaGetLastError", __FILE__, __LINE__);


    // diagonalize A
    int *d_info = NULL;
    double *d_Omega = NULL;
    check_Cuda_Errors(cudaMalloc((void**)&d_info, sizeof(int)),
        "cudaMalloc", __FILE__, __LINE__);
    check_Cuda_Errors(cudaMalloc((void**)&d_Omega, nS * sizeof(double)),
        "cudaMalloc", __FILE__, __LINE__);

    diag_dn_dsyevd(nS, d_info, d_Omega, d_A);
    check_Cuda_Errors(cudaGetLastError(), "cudaGetLastError", __FILE__, __LINE__);

    int info_gpu = 0;
    check_Cuda_Errors(cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost),
        "cudaMemcpy", __FILE__, __LINE__);
    if (info_gpu != 0) {
        printf("Error: diag_dn_dsyevd returned error code %d\n", info_gpu);
        exit(EXIT_FAILURE);
    }


    check_Cuda_Errors(cudaMemcpy(h_X, d_A, nS * nS * sizeof(double), cudaMemcpyDeviceToHost), 
        "cudaMemcpy", __FILE__, __LINE__);

    check_Cuda_Errors(cudaMemcpy(h_Omega, d_Omega, nS * sizeof(double), cudaMemcpyDeviceToHost), 
        "cudaMemcpy", __FILE__, __LINE__);

    check_Cuda_Errors(cudaFree(d_info), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_eps), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_ERI), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_A), "cudaFree", __FILE__, __LINE__);
    check_Cuda_Errors(cudaFree(d_Omega), "cudaFree", __FILE__, __LINE__);


}

