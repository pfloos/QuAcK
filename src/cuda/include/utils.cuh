#ifndef UTILS
#define UTILS

extern "C" void check_Cuda_Errors(cudaError_t err, const char* msg, const char* file, int line);

extern "C" void check_Cublas_Errors(cublasStatus_t status, const char* msg, const char* file, int line);


#endif
