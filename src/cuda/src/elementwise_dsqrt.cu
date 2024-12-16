#include <stdio.h>
#include <math.h>


__global__ void elementwise_dsqrt_kernel(int nS, double *A, double *A_Sq) {


    int i;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    while(i < nS) {

        if(A[i] > 0.0) {

            A_Sq[i] = sqrt(A[i]);

        } else {

            A_Sq[i] = sqrt(-A[i]);

        }

        i += blockDim.x * gridDim.x;
    } // i

}





extern "C" void elementwise_dsqrt(int nS, double *A, double *A_Sq) {

    int sBlocks = 32;
    int nBlocks = (nS + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, 1, 1);
    dim3 dimBlock(sBlocks, 1, 1);

    printf("lunching elementwise_dsqrt_kernel with %d blocks and %d threads/block\n",
        nBlocks, sBlocks);


    elementwise_dsqrt_kernel<<<dimGrid, dimBlock>>>(nS, A, A_Sq);

}




