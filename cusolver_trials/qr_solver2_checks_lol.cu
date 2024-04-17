#include <iostream>
#include <iomanip>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cublas_v2.h>

#define BLOCK_SIZE 32

// Error checking macro for CUDA calls
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// Error handling for cuSolver calls
inline void cusolveSafeCall(cusolverStatus_t status) {
    if(status != CUSOLVER_STATUS_SUCCESS) {
        std::cerr << "cuSolver API failed with status " << status << std::endl;
        exit(1);
    }
}

// Error handling for cuBLAS calls
inline void cublasSafeCall(cublasStatus_t status) {
    if(status != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "cuBLAS API failed with status " << status << std::endl;
        exit(1);
    }
}

__global__ void printMatrix(const double *A, int numRows, int numCols) {
    int row = threadIdx.x + blockIdx.x * blockDim.x;
    int col = threadIdx.y + blockIdx.y * blockDim.y;

    if(row < numRows && col < numCols) {
        printf("A[%d, %d] = %f\n", row, col, A[row + col * numRows]);
    }
}

void solver( double* A, double* b, double* x, int m, int n){
   // cuSOLVER and CUBLAS initialization
    cusolverDnHandle_t solver_handle;
    cublasHandle_t cublas_handle;
    cusolveSafeCall(cusolverDnCreate(&solver_handle));
    cublasSafeCall(cublasCreate(&cublas_handle));

    double *d_A, *d_TAU, *d_b, *work;
    int *devInfo, work_size = 0;

    gpuErrchk(cudaMalloc((void**)&d_A, m * n * sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&d_b, m * sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&d_TAU, n * sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&devInfo, sizeof(int)));

    gpuErrchk(cudaMemcpy(d_A, A, m * n * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_b, b, m * sizeof(double), cudaMemcpyHostToDevice));

    cusolveSafeCall(cusolverDnDgeqrf_bufferSize(solver_handle, m, n, d_A, m, &work_size));
    gpuErrchk(cudaMalloc((void**)&work, work_size * sizeof(double)));

    // QR decomposition 
    cusolveSafeCall(cusolverDnDgeqrf(solver_handle, m, n, d_A, m, d_TAU, work, work_size, devInfo));

    // Extracting the R matrix and compute Q^T*b
    double *d_R;
    gpuErrchk(cudaMalloc(&d_R, n * n * sizeof(double))); // R is NxN (upper triangular)
    dim3 Grid(BLOCK_SIZE, BLOCK_SIZE);
    dim3 Block((n + BLOCK_SIZE - 1) / BLOCK_SIZE, (n + BLOCK_SIZE - 1) / BLOCK_SIZE);

    printMatrix<<<Block, Grid>>>(d_A, m, n);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    // Computing Q^T*b (storing the result in d_b)
    cusolveSafeCall(cusolverDnDormqr(solver_handle, CUBLAS_SIDE_LEFT, CUBLAS_OP_T, m, 1, n, d_A, m, d_TAU, d_b, m, work, work_size, devInfo));

    // Solving R*x = Q^T*b for x
    const double alpha = 1.0;
    cublasSafeCall(cublasDtrsm(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, 1, &alpha, d_A, m, d_b, m));

    // Copying the d_b solution into h_x for printing
    double h_x[n];
    gpuErrchk(cudaMemcpy(h_x, d_b, n * sizeof(double), cudaMemcpyDeviceToHost));

    // Printing to terminal x[i] solution
    std::cout << "Solution x: \n";
    for(int i = 0; i < n; i++) {
        std::cout << "x[" << i << "] = " << h_x[i] << std::endl;
    }

    // DESTROY DEMOLISH DEVOUR 
    cudaFree(d_A);
    cudaFree(d_TAU);
    cudaFree(d_b);
    cudaFree(work);
    cudaFree(devInfo);
    cudaFree(d_R);

    cusolverDnDestroy(solver_handle);
    cublasDestroy(cublas_handle);

    // Cleanup
    cudaFree(d_A);
    cudaFree(d_tau);
    cudaFree(d_b);
    cudaFree(devInfo);
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    cublasDestroy(cublasH);
}

int main()
{
    int m = 3, n = 3;
    double A[m*n] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
    double b[m] = {1, 2, 3};
    double x[m];
    solver( A, b, x, m, n);
    for(int i =0;i<m;++i)
    {
        printf("%f\n",x[i]);
    }
}
/* double h_cuSolverMatrix[] = {0.0860585,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,1,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,1,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6.27048,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,20.2695,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5.29894,-4.5994,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,16.1788,15.3409,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,-1.18607,-2.02373,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,4.99719,5.16531,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,-1.76187,-2.00586,-6.97782,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,5.34622,5.11969,28.9878,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.08756,-1.77754,0,-6.97098,-2.955e-309,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,4.4891,5.42773,0,21.989,4.24576,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-2.81706e-309,-6.93031,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,4.04736,14.0862,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.6078e-309,-3.12449e-312,5.69561e-309,-1.16359e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,4.87286,0,2.82407e-309,15.787,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16193e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,6.02897,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16255e-308,0,5.92739,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,1.83216,0,2.82407e-309,9.65358,0,18.3623,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,7.38315,-6.87954,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,15.2761,2.56435,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,-2.09185,-1.16091e-308,0,11.5037,-6.86825,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,4.38058,0,0,4.33066,2.53794,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,-1.63282,-1.16091e-308,0,11.5037,-4.89132,0,-2.61345,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,3.39832,0,0,4.33066,-2.08962,0,2.3646,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,-3.31593,-1.16091e-308,0,11.5037,-4.89132,0,-2.62589,-1.17793,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,6.53021,0,0,4.33066,-2.08962,0,2.37586,2.39834,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,9.99991,-4.89132,0,-5.69091e-310,-1.19385,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,7.52928,-2.08962,0,-6.28982e-310,2.43075,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,3.59953e-309,-3.12449e-312,5.69561e-309,-1.16091e-308,0,11.5037,-4.89132,0,-5.69091e-310,-1.59089e-321,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,2.82407e-309,0,0,4.33066,-2.08962,0,-6.28982e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,-6.27048,-5.29894,3.30357e-312,3.30357e-312,-1.08756,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.99431e-312,1,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,-4.5994,-1.18607,-1.76187,-1.77754,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-2.02373,-2.00586,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,-6.97782,-6.97098,3.59953e-309,3.6078e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,3.59953e-309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-2.955e-309,-2.81706e-309,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,-3.12449e-312,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,-6.93031,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,-2.09185,-1.63282,-3.31593,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,5.69561e-309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16359e-308,-1.16193e-308,-1.16255e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,-1.16091e-308,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,5.92739,7.38315,11.5037,11.5037,11.5037,9.99991,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,11.5037,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-6.87954,-6.86825,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,-4.89132,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-2.61345,-2.62589,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,-5.69091e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.17793,-1.19385,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,-1.59089e-321,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,20.2695,16.1788,8.00344e-313,8.00344e-313,4.4891,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,15.3409,4.99719,5.34622,5.42773,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.16531,5.11969,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,28.9878,21.989,5.13e-312,4.87286,5.13e-312,1.83216,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.24576,4.04736,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,14.0862,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,4.38058,3.39832,6.53021,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,2.82407e-309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15.787,6.02897,9.65358,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,18.3623,15.2761,4.33066,4.33066,4.33066,7.52928,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,4.33066,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,2.56435,2.53794,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,-2.08962,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,2.3646,2.37586,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,-6.28982e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,2.39834,2.43075,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,6.92437e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double h_cuSolverVector[] = {40,40,40,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,225.753,13.4936,-96.244,-54.9104,-7.6,5.38456,-2.84615e-311,24,-35.3885,-11.0896,-3.5,-8.85719,-14.7368,-14.9,25.3857,14.4312,-8.28304,29.5385,-1.66932,-1.52778,10.9658,-1.8,13.2174,-0.908676,-5,2.49465,2.51826,-4.6703e-310,0,0,0,0,0,0,0,0,0,0,0,0};*/
