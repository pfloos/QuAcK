#include <cuda_runtime.h>
#include <cuda.h>
#include <stdio.h>
#include <cublas_v2.h>


extern "C" void check_Cuda_Errors(cudaError_t err, const char* msg, const char* file, int line) {
    if (err != cudaSuccess) {
        printf("CUDA Error in %s at line %d\n", file, line);
        printf("%s - %s\n", msg, cudaGetErrorString(err));
        exit(0);
    }
}


const char* cublasGetErrorString(cublasStatus_t status) {
    switch (status) {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
        case CUBLAS_STATUS_NOT_SUPPORTED:
            return "CUBLAS_STATUS_NOT_SUPPORTED";
        case CUBLAS_STATUS_LICENSE_ERROR:
            return "CUBLAS_STATUS_LICENSE_ERROR";
    }
    return "UNKNOWN CUBLAS ERROR";
}

extern "C" void check_Cublas_Errors(cublasStatus_t status, const char* msg, const char* file, int line) {

    const char* err = cublasGetErrorString(status);

    if (err != "CUBLAS_STATUS_SUCCESS") {
        printf("CUBLAS Error in %s at line %d\n", file, line);
        printf("%s - %s\n", msg, err);
        exit(0);
    }
}


