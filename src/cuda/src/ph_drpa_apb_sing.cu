#include <stdio.h>

__global__ void ph_dRPA_ApB_sing_kernel(int nO, int nV, int nBas, int nS, double *eps, double *ERI, double *ApB) {


    int i, j, a, b;
    int aa, bb;
    int nVS;
    int nBas2, nBas3;
    int i_A0, i_A1, i_A2;
    int i_I0, i_I1, i_I2;

    bool a_eq_b;

    nVS = nV * nS;

    nBas2 = nBas * nBas;
    nBas3 = nBas2 * nBas;

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

            i = 0;
            while(i < nO) {

                i_A2 = i_A1 + i * nVS;
                i_I2 = i_I1 + i;
 
                j = 0;
                while(j < nO) {

                    ApB[i_A2 + j * nV] = 2.0 * (ERI[i_I2 + j * nBas3] + ERI[i_I2 + j * nBas]);
                    if(a_eq_b && (i==j)) {
                        ApB[i_A2 + j * nV] += eps[a] - eps[i];
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




