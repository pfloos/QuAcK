#include <stdio.h>

__global__ void ph_dRPA_B_sing_kernel(int nO, int nV, int nBas, int nS, double *ERI, double *B) {


    int i, j, a, b;
    int aa, bb;
    int nVS;
    int nBas2, nBas3;
    int i_B0, i_B1, i_B2;
    int i_I0, i_I1, i_I2;

    nVS = nV * nS;

    nBas2 = nBas * nBas;
    nBas3 = nBas2 * nBas;

    aa = blockIdx.x * blockDim.x + threadIdx.x;
    bb = blockIdx.y * blockDim.y + threadIdx.y;

    while(aa < nV) {
        a = aa + nO;

        i_B0 = aa * nS;
        i_I0 = a * nBas2;

        while(bb < nV) {
            b = bb + nO;

            i_B1 = i_B0 + bb;
            i_I1 = i_I0 + b * nBas3;

            i = 0;
            while(i < nO) {

                i_B2 = i_B1 + i * nVS;
                i_I2 = i_I1 + i;
 
                j = 0;
                while(j < nO) {

                    B[i_B2 + j * nV] = 2.0 * ERI[i_I2 + j * nBas];

                    j ++;
	        } // j

                i ++;
            } // i

            bb += blockDim.y * gridDim.y;
        } // bb

        aa += blockDim.x * gridDim.x;
    } // aa

}





extern "C" void ph_dRPA_B_sing(int nO, int nV, int nBas, int nS, double *ERI, double *B) {


    int sBlocks = 32;
    int nBlocks = (nV + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, nBlocks, 1);
    dim3 dimBlock(sBlocks, sBlocks, 1);


    printf("lunching ph_dRPA_B_sing_kernel with %dx%d blocks and %dx%d threads/block\n",
        nBlocks, nBlocks, sBlocks, sBlocks);


    ph_dRPA_B_sing_kernel<<<dimGrid, dimBlock>>>(nO, nV, nBas, nS, ERI, B);

}




