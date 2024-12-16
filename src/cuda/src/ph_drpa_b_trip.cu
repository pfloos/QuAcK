#include <stdio.h>

__global__ void ph_dRPA_B_trip_kernel(int nO, int nV, int nBas, int nS, double *B) {


    int i, j;
    int aa, bb;
    int nVS;
    int i_B0, i_B1, i_B2;

    nVS = nV * nS;

    aa = blockIdx.x * blockDim.x + threadIdx.x;
    bb = blockIdx.y * blockDim.y + threadIdx.y;

    while(aa < nV) {

        i_B0 = aa * nS;

        while(bb < nV) {

            i_B1 = i_B0 + bb;

            i = 0;
            while(i < nO) {

                i_B2 = i_B1 + i * nVS;
 
                j = 0;
                while(j < nO) {

                    B[i_B2 + j * nV] = 0.0;

                    j ++;
	        } // j

                i ++;
            } // i

            bb += blockDim.y * gridDim.y;
        } // bb

        aa += blockDim.x * gridDim.x;
    } // aa

}





extern "C" void ph_dRPA_B_trip(int nO, int nV, int nBas, int nS, double *B) {


    int sBlocks = 32;
    int nBlocks = (nV + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, nBlocks, 1);
    dim3 dimBlock(sBlocks, sBlocks, 1);


    printf("lunching ph_dRPA_B_trip_kernel with %dx%d blocks and %dx%d threads/block\n",
        nBlocks, nBlocks, sBlocks, sBlocks);


    ph_dRPA_B_trip_kernel<<<dimGrid, dimBlock>>>(nO, nV, nBas, nS, B);

}




