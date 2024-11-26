#ifndef PH_DRPA

#define PH_DRPA

extern void check_Cuda_Errors(cudaError_t err, const char * msg, const char * file, int line);
extern void check_Cublas_Errors(cublasStatus_t status, const char * msg, const char * file, int line);

extern void phLR_dRPA_A_sing(int nO, int nBas, double *eps, double *ERI, double *A);

#endif
