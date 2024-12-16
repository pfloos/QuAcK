#include <stdio.h>


__global__ void ph_dRPA_ApB_sing_kernel(int nO, int nV, int nBas, int nS, 
                                        double *eps, double *ERI, double *ApB) {


    long i, j, a, b;
    long aa, bb;

    int i_A0, i_A1, i_A2, i_A3;
    int i_I0, i_I1, i_I2;
    int i_J1, i_J2;

    int nVS;
    int nBas2;

    long long i_I3, i_J3;
    long long nBas3;

    bool a_eq_b;

    nVS = nV * nS;

    nBas2 = nBas * nBas;
    nBas3 = (long long) nBas2 * (long long) nBas;

    aa = blockIdx.x * blockDim.x + threadIdx.x;
    bb = blockIdx.y * blockDim.y + threadIdx.y;

    while(aa < nV) {
        a = aa + nO;

        i_A0 = aa * nS;
        i_I0 = a * nBas2;

        while(bb < nV) {
            b = bb + nO;

            a_eq_b = a == b;

            i_A1 = i_A0 + bb;
            i_I1 = i_I0 + b * nBas;
            i_J1 = a + b * nBas;

            i = 0;
            while(i < nO) {

                i_A2 = i_A1 + i * nVS;
                i_I2 = i_I1 + i;
                i_J2 = i_J1 + i * nBas2;
 
                j = 0;
                while(j < nO) {

                    i_A3 = i_A2 + j * nV;
                    i_I3 = i_I2 + (long long) j * nBas3;
                    i_J3 = i_J2 + (long long) j * nBas3;

                    ApB[i_A3] = 2.0 * (ERI[i_I3] + ERI[i_J3]);
                    if(a_eq_b && (i==j)) {
                        ApB[i_A3] += eps[a] - eps[i];
                    }

                    j ++;
	        } // j

                i ++;
            } // i

            bb += blockDim.y * gridDim.y;
        } // bb

        aa += blockDim.x * gridDim.x;
    } // aa

}





extern "C" void ph_dRPA_ApB_sing(int nO, int nV, int nBas, int nS, double *eps, double *ERI, double *ApB) {


    int sBlocks = 32;
    int nBlocks = (nV + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, nBlocks, 1);
    dim3 dimBlock(sBlocks, sBlocks, 1);

    printf("lunching ph_dRPA_ApB_sing_kernel with %dx%d blocks and %dx%d threads/block\n",
        nBlocks, nBlocks, sBlocks, sBlocks);


    ph_dRPA_ApB_sing_kernel<<<dimGrid, dimBlock>>>(nO, nV, nBas, nS, eps, ERI, ApB);

}




