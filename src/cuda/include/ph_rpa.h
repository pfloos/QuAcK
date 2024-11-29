#ifndef PH_RPA

#define PH_RPA

extern void ph_dRPA_A_sing(int nO, int nV, int nBas, int nS, double *eps, double *ERI, double *A);
extern void ph_dRPA_B_sing(int nO, int nV, int nBas, int nS, double *ERI, double *B);

extern void ph_dRPA_ApB_sing(int nO, int nV, int nBas, int nS, double *eps, double *ERI, double *ApB);
extern void ph_dRPA_AmB_sing(int nO, int nV, int nBas, int nS, double *eps, double *ERI, double *AmB);

#endif
