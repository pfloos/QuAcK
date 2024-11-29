#include <stdio.h>

__global__ void ph_dRPA_B_sing_kernel(int nO, int nV, int nBas, int nS, double *ERI, double *B) {


    int i, j, a, b;
    int aa, bb;

    long long nVS;
    long long nBas2, nBas3;
    long long i_B0, i_B1, i_B2, i_B3;
    long long i_I0, i_I1, i_I2, i_I3;


    nVS = (long long) nV * (long long) nS;

    nBas2 = (long long) nBas * (long long) nBas;
    nBas3 = nBas2 * (long long) nBas;

    aa = blockIdx.x * blockDim.x + threadIdx.x;
    bb = blockIdx.y * blockDim.y + threadIdx.y;

    while(aa < nV) {
        a = aa + nO;

        i_B0 = (long long) aa * (long long) nS;
        i_I0 = (long long) a * nBas2;

        while(bb < nV) {
            b = bb + nO;

            i_B1 = i_B0 + (long long) bb;
            i_I1 = i_I0 + (long long) b * nBas3;

            i = 0;
            while(i < nO) {

                i_B2 = i_B1 + (long long) i * nVS;
                i_I2 = i_I1 + (long long) i;
 
                j = 0;
                while(j < nO) {

                    i_B3 = i_B2 + (long long) j * (long long) nV;
                    i_I3 = i_I2 + (long long) j * (long long) nBas;

                    B[i_B3] = 2.0 * ERI[i_I3];

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




