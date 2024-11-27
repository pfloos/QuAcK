#ifndef UTILS

#define UTILS

extern void check_Cuda_Errors(cudaError_t err, const char *msg, const char *file, int line);
extern void check_Cublas_Errors(cublasStatus_t status, const char *msg, const char *file, int line);
extern void check_Cusolver_Errors(cusolverStatus_t status, const char *msg, const char *file, int line);

#endif
