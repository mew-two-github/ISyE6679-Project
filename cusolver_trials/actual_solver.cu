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

int main() {
//                    /\_[]_/\
//               |] _||_ [|
//        ___     \/ || \/
//       /___\       ||
//      (|0 0|)      ||
//    __/{\U/}\_ ___/vvv
//   / \  {~}   / _|_P|
//   | /\  ~   /_/   []
//   |_| (____)        
//   \_]/______\        -Defending the gates to the most flawlessly written cuSolver code-
//      _\_||_/_           
//     (_,_||_,_)

    const int Nrows = 3;
    const int Ncols = 3;

    //ASSUMPTION: Nrows >= Ncols 
    double h_A[Nrows * Ncols] = {1, 4, 7, 2, 5, 8, 3, 6, 10};  // Column-major storage
    double h_b[Nrows] = {1, 2, 3};

    // cuSOLVER and CUBLAS initialization
    cusolverDnHandle_t solver_handle;
    cublasHandle_t cublas_handle;
    cusolveSafeCall(cusolverDnCreate(&solver_handle));
    cublasSafeCall(cublasCreate(&cublas_handle));

    double *d_A, *d_TAU, *d_b, *work;
    int *devInfo, work_size = 0;

    gpuErrchk(cudaMalloc((void**)&d_A, Nrows * Ncols * sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&d_b, Nrows * sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&d_TAU, Ncols * sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&devInfo, sizeof(int)));

    gpuErrchk(cudaMemcpy(d_A, h_A, Nrows * Ncols * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_b, h_b, Nrows * sizeof(double), cudaMemcpyHostToDevice));

    cusolveSafeCall(cusolverDnDgeqrf_bufferSize(solver_handle, Nrows, Ncols, d_A, Nrows, &work_size));
    gpuErrchk(cudaMalloc((void**)&work, work_size * sizeof(double)));

    // QR decomposition 
    cusolveSafeCall(cusolverDnDgeqrf(solver_handle, Nrows, Ncols, d_A, Nrows, d_TAU, work, work_size, devInfo));

    // Extracting the R matrix and compute Q^T*b
    double *d_R;
    gpuErrchk(cudaMalloc(&d_R, Ncols * Ncols * sizeof(double))); // R is NxN (upper triangular)
    dim3 Grid(BLOCK_SIZE, BLOCK_SIZE);
    dim3 Block((Ncols + BLOCK_SIZE - 1) / BLOCK_SIZE, (Ncols + BLOCK_SIZE - 1) / BLOCK_SIZE);

    printMatrix<<<Block, Grid>>>(d_A, Nrows, Ncols);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    // Computing Q^T*b (storing the result in d_b)
    cusolveSafeCall(cusolverDnDormqr(solver_handle, CUBLAS_SIDE_LEFT, CUBLAS_OP_T, Nrows, 1, Ncols, d_A, Nrows, d_TAU, d_b, Nrows, work, work_size, devInfo));

    // Solving R*x = Q^T*b for x
    const double alpha = 1.0;
    cublasSafeCall(cublasDtrsm(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, Ncols, 1, &alpha, d_A, Nrows, d_b, Nrows));

    // Copying the d_b solution into h_x for printing
    double h_x[Ncols];
    gpuErrchk(cudaMemcpy(h_x, d_b, Ncols * sizeof(double), cudaMemcpyDeviceToHost));

    // Printing to terminal x[i] solution
    std::cout << "Solution x: \n";
    for(int i = 0; i < Ncols; i++) {
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

    return 0;
}
