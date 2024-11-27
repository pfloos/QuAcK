#ifndef PH_DRPA

#define PH_DRPA

extern void check_Cuda_Errors(cudaError_t err, const char * msg, const char * file, int line);
extern void check_Cublas_Errors(cublasStatus_t status, const char * msg, const char * file, int line);

extern void ph_dRPA_A_sing(int nO, int nV, int nBas, int nS, double *eps, double *ERI, double *A);

#endif
