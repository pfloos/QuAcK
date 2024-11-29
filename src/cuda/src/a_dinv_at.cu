#include <stdio.h>


__global__ void A_Dinv_At_kernel(int n, double *A, double *D, double *R) {


    int i, j;
    int k;
    int in, ij;
    int kn;

    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;

    while(i < n) {

        in = i * n;

        while(j < n) {

            ij = in + j;

            R[ij] = 0.0;
            k = 0;
            while(k < n) {

                kn = k * n;
                R[ij] += D[k] * U[i + kn] * U[j + kn] / (D[k] + 1e-12);

                k ++;
            } // k

            j += blockDim.y * gridDim.y;
        } // j

        i += blockDim.x * gridDim.x;
    } // i

}





extern "C" void A_Dinv_At(int n, double *A, double *D, double *R) {


    int sBlocks = 32;
    int nBlocks = (n + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, nBlocks, 1);
    dim3 dimBlock(sBlocks, sBlocks, 1);

    printf("lunching A_Dinv_At_kernel with %dx%d blocks and %dx%d threads/block\n",
        nBlocks, nBlocks, sBlocks, sBlocks);


    A_Dinv_At_kernel<<<dimGrid, dimBlock>>>(n, A, D, R);

}




