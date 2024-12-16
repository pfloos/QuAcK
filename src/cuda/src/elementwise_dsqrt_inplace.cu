#include <stdio.h>
#include <math.h>


__global__ void elementwise_dsqrt_inplace_kernel(int n, double *A) {


    int i;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    while(i < n) {

        if(A[i] > 0.0) {

            A[i] = sqrt(A[i]);

        } else {

            A[i] = sqrt(-A[i]);

        }

        i += blockDim.x * gridDim.x;
    } // i

}





extern "C" void elementwise_dsqrt_inplace(int n, double *A) {

    int sBlocks = 32;
    int nBlocks = (n + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, 1, 1);
    dim3 dimBlock(sBlocks, 1, 1);

    printf("lunching elementwise_dsqrt_inplace_kernel with %d blocks and %d threads/block\n",
        nBlocks, sBlocks);


    elementwise_dsqrt_inplace_kernel<<<dimGrid, dimBlock>>>(n, A);

}




