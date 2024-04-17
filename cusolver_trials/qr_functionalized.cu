#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cublas_v2.h>

// Error handling macros used by CUDA runtime methods
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// Error handling for cuSolver operations
inline void cusolverSafeCall(cusolverStatus_t status) {
    if(status != CUSOLVER_STATUS_SUCCESS) {
        std::cerr << "cuSolver API failed with status " << status << std::endl;
        exit(1);
    }
}

// Error handling for cuBLAS operations
inline void cublasSafeCall(cublasStatus_t status) {
    if(status != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "cuBLAS API failed with status " << status << std::endl;
        exit(1);
    }
}

void solver( double* A, double* b, double* x, int Nrows, int Ncols) {
    cusolverDnHandle_t solver_handle;
    cublasHandle_t cublas_handle;
    cusolverSafeCall(cusolverDnCreate(&solver_handle));
    cublasSafeCall(cublasCreate(&cublas_handle));

    double *d_A, *d_TAU, *d_b, *work;
    int *devInfo, work_size = 0;

    gpuErrchk(cudaMalloc((void**)&d_A, Nrows * Ncols * sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&d_b, Nrows * sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&d_TAU, Ncols * sizeof(double)));
    gpuErrchk(cudaMalloc((void**)&devInfo, sizeof(int)));

    gpuErrchk(cudaMemcpy(d_A, h_A.data(), Nrows * Ncols * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_b, h_b.data(), Nrows * sizeof(double), cudaMemcpyHostToDevice));

    cusolverSafeCall(cusolverDnDgeqrf_bufferSize(solver_handle, Nrows, Ncols, d_A, Nrows, &work_size));
    gpuErrchk(cudaMalloc((void**)&work, work_size * sizeof(double)));

    cusolverSafeCall(cusolverDnDgeqrf(solver_handle, Nrows, Ncols, d_A, Nrows, d_TAU, work, work_size, devInfo));
    cusolverSafeCall(cusolverDnDormqr(solver_handle, CUBLAS_SIDE_LEFT, CUBLAS_OP_T, Nrows, 1, Ncols, d_A, Nrows, d_TAU, d_b, Nrows, work, work_size, devInfo));

    const double alpha = 1.0;
    cublasSafeCall(cublasDtrsm(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, Ncols, 1, &alpha, d_A, Nrows, d_b, Nrows));

    std::vector<double> h_x(Ncols);
    gpuErrchk(cudaMemcpy(h_x.data(), d_b, Ncols * sizeof(double), cudaMemcpyDeviceToHost));

    cudaFree(d_A);
    cudaFree(d_TAU);
    cudaFree(d_b);
    cudaFree(work);
    cudaFree(devInfo);

    cusolverDnDestroy(solver_handle);
    cublasDestroy(cublas_handle);
}

int main() {
    int m = 3, n = 3;
    double A[m*n] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
    double b[m] = {1, 2, 3};
    double x[m];

    solveQR(A, b, m, n);

    std::cout << "Solution x: \n";
    for (int i = 0; i < result.size(); i++) {
        std::cout << "x[" << i << "] = " << result[i] << std::endl;
    }

    return 0;
}
