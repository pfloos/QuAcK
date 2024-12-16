#include <stdio.h>

__global__ void ph_dRPA_A_trip_kernel(int nO, int nV, int nBas, int nS, double *eps, double *A) {


    int i, j, a, b;
    int aa, bb;
    int nVS;
    int i_A0, i_A1, i_A2;

    nVS = nV * nS;

    aa = blockIdx.x * blockDim.x + threadIdx.x;
    bb = blockIdx.y * blockDim.y + threadIdx.y;

    while(aa < nV) {
        a = aa + nO;

        i_A0 = aa * nS;

        while(bb < nV) {
            b = bb + nO;

            i_A1 = i_A0 + bb;

            i = 0;
            while(i < nO) {

                i_A2 = i_A1 + i * nVS;
 
                j = 0;
                while(j < nO) {

                    A[i_A2 + j * nV] = 0.0;
                    if((a==b) && (i==j)) {
                        A[i_A2 + j * nV] += eps[a] - eps[i];
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





extern "C" void ph_dRPA_A_trip(int nO, int nV, int nBas, int nS, double *eps, double *A) {


    int sBlocks = 32;
    int nBlocks = (nV + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, nBlocks, 1);
    dim3 dimBlock(sBlocks, sBlocks, 1);


    printf("lunching ph_dRPA_A_trip_kernel with %dx%d blocks and %dx%d threads/block\n",
        nBlocks, nBlocks, sBlocks, sBlocks);


    ph_dRPA_A_trip_kernel<<<dimGrid, dimBlock>>>(nO, nV, nBas, nS, eps, A);

}




