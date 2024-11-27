#include <cuda_runtime.h>
#include <cuda.h>
#include <stdio.h>
#include <cublas_v2.h>
#include <cstring>
#include <cusolverDn.h>



extern "C" void check_Cuda_Errors(cudaError_t err, const char* msg, const char* file, int line) {
    if (err != cudaSuccess) {
        printf("CUDA Error in %s at line %d\n", file, line);
        printf("%s - %s\n", msg, cudaGetErrorString(err));
        exit(0);
    }
}


const char* cublas_Get_Error_String(cublasStatus_t status) {
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

    const char* err = cublas_Get_Error_String(status);

    if (strcmp(err, "CUBLAS_STATUS_SUCCESS") != 0) {
        printf("CUBLAS Error in %s at line %d\n", file, line);
        printf("%s - %s\n", msg, err);
        exit(0);
    }
}



const char* cusolver_Get_Error_String(cusolverStatus_t status) {
    switch (status) {
        case CUSOLVER_STATUS_SUCCESS:
            return "CUSOLVER_STATUS_SUCCESS";
        case CUSOLVER_STATUS_NOT_INITIALIZED:
            return "CUSOLVER_STATUS_NOT_INITIALIZED";
        case CUSOLVER_STATUS_ALLOC_FAILED:
            return "CUSOLVER_STATUS_ALLOC_FAILED";
        case CUSOLVER_STATUS_INVALID_VALUE:
            return "CUSOLVER_STATUS_INVALID_VALUE";
        case CUSOLVER_STATUS_ARCH_MISMATCH:
            return "CUSOLVER_STATUS_ARCH_MISMATCH";
        case CUSOLVER_STATUS_MAPPING_ERROR:
            return "CUSOLVER_STATUS_MAPPING_ERROR";
        case CUSOLVER_STATUS_EXECUTION_FAILED:
            return "CUSOLVER_STATUS_EXECUTION_FAILED";
        case CUSOLVER_STATUS_INTERNAL_ERROR:
            return "CUSOLVER_STATUS_INTERNAL_ERROR";
        case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
        case CUSOLVER_STATUS_NOT_SUPPORTED:
            return "CUSOLVER_STATUS_NOT_SUPPORTED";
        case CUSOLVER_STATUS_ZERO_PIVOT:
            return "CUSOLVER_STATUS_ZERO_PIVOT";
        case CUSOLVER_STATUS_INVALID_LICENSE:
            return "CUSOLVER_STATUS_INVALID_LICENSE";
        case CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED:
            return "CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED";
        case CUSOLVER_STATUS_IRS_PARAMS_INVALID:
            return "CUSOLVER_STATUS_IRS_PARAMS_INVALID";
        case CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC:
            return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC";
        case CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE:
            return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE";
        case CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER:
            return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER";
        case CUSOLVER_STATUS_IRS_INTERNAL_ERROR:
            return "CUSOLVER_STATUS_IRS_INTERNAL_ERROR";
        case CUSOLVER_STATUS_IRS_NOT_SUPPORTED:
            return "CUSOLVER_STATUS_IRS_NOT_SUPPORTED";
        case CUSOLVER_STATUS_IRS_OUT_OF_RANGE:
            return "CUSOLVER_STATUS_IRS_OUT_OF_RANGE";
        case CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES:
            return "CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES";
        case CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED:
            return "CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED";
        case CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED:
            return "CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED";
        case CUSOLVER_STATUS_IRS_MATRIX_SINGULAR:
            return "CUSOLVER_STATUS_IRS_MATRIX_SINGULAR";
        case CUSOLVER_STATUS_INVALID_WORKSPACE:
            return "CUSOLVER_STATUS_INVALID_WORKSPACE";
        default:
            return "UNKNOWN CUSOLVER ERROR";
    }
}

extern "C" void check_Cusolver_Errors(cusolverStatus_t status, const char* msg, const char* file, int line) {

    const char* err = cusolver_Get_Error_String(status);

    if (status != CUSOLVER_STATUS_SUCCESS) {
        printf("CUSOLVER Error in %s at line %d\n", file, line);
        printf("%s - %s\n", msg, err);
        exit(EXIT_FAILURE);
    }
}

