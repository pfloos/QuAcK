#include <stdio.h>


__global__ void A_minus_twoB_in_B_kernel(int n, double *A, double *B) {


    int i, j;
    int in, ji;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;

    while(i < n) {

        in = i * n;

        while(j < n) {

            ji = in + j;

            B[ji] = A[ji] - 2.0 * B[ji];

            j += blockDim.y * gridDim.y;
        } // j

        i += blockDim.x * gridDim.x;
    } // i

}




extern "C" void A_minus_twoB_in_B(int n, double *A, double *B) {


    int sBlocks = 32;
    int nBlocks = (n + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, nBlocks, 1);
    dim3 dimBlock(sBlocks, sBlocks, 1);

    printf("lunching A_minus_twoB_in_B_kernel with %dx%d blocks and %dx%d threads/block\n",
        nBlocks, nBlocks, sBlocks, sBlocks);


    A_minus_twoB_in_B_kernel<<<dimGrid, dimBlock>>>(n, A, B);

}



