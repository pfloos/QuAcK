#include <stdio.h>


__global__ void A_D_in_B_kernel(int n, double *A, double *D, double *B) {


    int i, j;
    int in, ji;

    double tmp;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;

    while(i < n) {

        in = i * n;

        tmp = D[i];

        while(j < n) {

            ji = in + j;

            B[ji] = A[ji] * tmp;

            j += blockDim.y * gridDim.y;
        } // j

        i += blockDim.x * gridDim.x;
    } // i

}





extern "C" void A_D_in_B(int n, double *A, double *D, double *B) {


    int sBlocks = 32;
    int nBlocks = (n + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, nBlocks, 1);
    dim3 dimBlock(sBlocks, sBlocks, 1);

    printf("lunching A_D_in_B_kernel with %dx%d blocks and %dx%d threads/block\n",
        nBlocks, nBlocks, sBlocks, sBlocks);


    A_D_in_B_kernel<<<dimGrid, dimBlock>>>(n, A, D, B);

}



