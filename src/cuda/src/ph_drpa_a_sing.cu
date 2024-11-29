#include <stdio.h>

__global__ void ph_dRPA_A_sing_kernel(int nO, int nV, int nBas, int nS, double *eps, double *ERI, double *A) {


    int i, j, a, b;
    int aa, bb;

    long long nVS;
    long long nBas2, nBas3;
    long long i_A0, i_A1, i_A2, i_A3;
    long long i_I0, i_I1, i_I2, i_I3;

    bool a_eq_b;

    nVS = (long long) nV * (long long) nS;

    nBas2 =  (long long) nBas *  (long long) nBas;
    nBas3 = nBas2 *  (long long) nBas;

    aa = blockIdx.x * blockDim.x + threadIdx.x;
    bb = blockIdx.y * blockDim.y + threadIdx.y;

    while(aa < nV) {
        a = aa + nO;

        i_A0 = (long long) aa * (long long) nS;
        i_I0 = (long long) a * nBas2;

        while(bb < nV) {
            b = bb + nO;

            a_eq_b = a == b;

            i_A1 = i_A0 + (long long) bb;
            i_I1 = i_I0 + (long long) b * (long long) nBas;

            i = 0;
            while(i < nO) {

                i_A2 = i_A1 + (long long) i * nVS;
                i_I2 = i_I1 + (long long) i;
 
                j = 0;
                while(j < nO) {

                    i_A3 = i_A2 + (long long) j * (long long) nV;
                    i_I3 = i_I2 + (long long) j * nBas3;

                    A[i_A3] = 2.0 * ERI[i_I3];
                    if(a_eq_b && (i==j)) {
                        A[i_A3] += eps[a] - eps[i];
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





extern "C" void ph_dRPA_A_sing(int nO, int nV, int nBas, int nS, double *eps, double *ERI, double *A) {


    int sBlocks = 32;
    int nBlocks = (nV + sBlocks - 1) / sBlocks;

    dim3 dimGrid(nBlocks, nBlocks, 1);
    dim3 dimBlock(sBlocks, sBlocks, 1);

    //dim3 dimGrid(nBlocks, 1, 1);
    //dim3 dimBlock(sBlocks, 1, 1);

    printf("lunching ph_dRPA_A_sing_kernel with %dx%d blocks and %dx%d threads/block\n",
        nBlocks, nBlocks, sBlocks, sBlocks);


    ph_dRPA_A_sing_kernel<<<dimGrid, dimBlock>>>(nO, nV, nBas, nS, eps, ERI, A);

}




