#ifndef MY_LINALG

#define MY_LINALG

extern void A_plus_B_in_A(int n, double *A, double *B);
extern void A_minus_twoB_in_B(int n, double *A, double *B);

extern void A_D_At(int n, double *A, double *D, double *R);
extern void A_Dinv_At(int n, double *A, double *D, double *R);

extern void A_D_inplace(int n, double *A, double *D);
extern void A_Dinv_inplace(int n, double *A, double *D);

extern void A_D_in_B(int n, double *A, double *D, double *B);
extern void A_Dinv_in_B(int n, double *A, double *D, double *B);

extern void elementwise_dsqrt(int nS, double *A, double *A_Sq);
extern void elementwise_dsqrt_inplace(int nS, double *A);

extern void diag_dn_dsyevd(int n, int *info, double *W, double *A);

#endif
