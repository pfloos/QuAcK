#include <stdio.h>

__global__ void phLR_dRPA_A_sing_kernel(int nO, int nBas, double *eps, double *ERI, double *A) {


    int i, j, a, b;
    int ia, jb, jb_off;

    int ij_off0, ij_off;

    int aa_max = nBas - nO;
    int ia_max = aa_max * nO;

    int nBas2 = nBas * nBas;
    int nBas3 = nBas2 * nBas;

    int aa = blockIdx.x * blockDim.x + threadIdx.x;
    int bb = blockIdx.y * blockDim.y + threadIdx.y;

    while(aa < aa_max) {
        a = aa + nO;

        ij_off0 = a * nBas2;

        while(bb < aa_max) {
            b = bb + nO;

            ij_off = ij_off0 + b * nBas;

            while(i < nO) {
                ia = i * aa_max + aa;
                jb_off = ia * ia_max;
 
                while(j < nO) {
                    jb = j * aa_max + bb;

                    A[jb + jb_off] = 2.0 * ERI[i + j * nBas3 + ij_off];
                    if(a==b && i==j) {
                        A[jb + jb_off] += eps[a] - eps[i];
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





extern "C" void phLR_dRPA_A_sing(int nO, int nBas, double *eps, double *ERI, double *A) {


    int size = nBas - nO;

    int sBlocks = 32;
    int nBlocks = (size + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, nBlocks, 1);
    dim3 dimBlock(sBlocks, sBlocks, 1);


    printf("lunching phLR_dRPA_A_sing_kernel with %dx%d blocks and %dx%d threads/block\n",
        nBlocks, nBlocks, sBlocks, sBlocks);


    phLR_dRPA_A_sing_kernel<<<dimGrid, dimBlock>>>(nO, nBas, eps, ERI, A);

}




