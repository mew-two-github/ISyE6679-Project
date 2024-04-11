#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cublas_v2.h>

// Error checking for CUDA calls
#define CHECK_CUDA(call) {\
    const cudaError_t error = call;\
    if (error != cudaSuccess) {\
        printf("Error: %s:%d, ", __FILE__, __LINE__);\
        printf("code:%d, reason: %s\n", error, cudaGetErrorString(error));\
        exit(1);\
    }\
}

// Error checking for cuSolver calls
#define CHECK_CUSOLVER(call, handle) {\
    const cusolverStatus_t error = call;\
    if (error != CUSOLVER_STATUS_SUCCESS) {\
        if (handle) cusolverDnDestroy(handle);\
        printf("Error: %s:%d, ", __FILE__, __LINE__);\
        printf("CUSOLVER error code:%d\n", error);\
        exit(1);\
    }\
}

// Error checking for cuBLAS calls
#define CHECK_CUBLAS(call, handle) {\
    const cublasStatus_t error = call;\
    if (error != CUBLAS_STATUS_SUCCESS) {\
        if (handle) cublasDestroy(handle);\
        printf("Error: %s:%d, ", __FILE__, __LINE__);\
        printf("CUBLAS error code:%d\n", error);\
        exit(1);\
    }\
}

int main() {
    cusolverDnHandle_t cusolverH = NULL;
    cublasHandle_t cublasH = NULL;
    const int m = 3; // Number of rows of A
    const int n = 3; // Number of columns of A
    double A[m*n] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
    double b[m] = {1, 2, 3}; 
    double *d_A = NULL, *d_tau = NULL, *d_b = NULL;
    int *devInfo = NULL;
    double *d_work = NULL;
    int lwork = 0;
    int info_gpu = 0;
    double alpha = 1.0;

    
    CHECK_CUDA(cudaSetDevice(0)); 
    CHECK_CUSOLVER(cusolverDnCreate(&cusolverH), cusolverH);
    CHECK_CUBLAS(cublasCreate(&cublasH), cublasH);

    // Allocate memory 
    CHECK_CUDA(cudaMalloc((void**)&d_A, sizeof(double) * m * n));
    CHECK_CUDA(cudaMalloc((void**)&d_tau, sizeof(double) * n));
    CHECK_CUDA(cudaMalloc((void**)&d_b, sizeof(double) * m));
    CHECK_CUDA(cudaMalloc((void**)&devInfo, sizeof(int)));

    // Copy host memory
    CHECK_CUDA(cudaMemcpy(d_A, A, sizeof(double) * m * n, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_b, b, sizeof(double) * m, cudaMemcpyHostToDevice));

    //space of geqrf
    CHECK_CUSOLVER(cusolverDnDgeqrf_bufferSize(cusolverH, m, n, d_A, m, &lwork), cusolverH);
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(double) * lwork));

    // QR factorization
    CHECK_CUSOLVER(cusolverDnDgeqrf(cusolverH, m, n, d_A, m, d_tau, d_work, lwork, devInfo), cusolverH);
    CHECK_CUDA(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));

    if (info_gpu != 0) {
        printf("QR factorization failed, info: %d\n", info_gpu);
    }

    // apply Q^T to b
    CHECK_CUSOLVER(cusolverDnDormqr(cusolverH, CUBLAS_SIDE_LEFT, CUBLAS_OP_T, m, 1, n, d_A, m, d_tau, d_b, m, d_work, lwork, devInfo), cusolverH);
    CHECK_CUDA(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));

    if (info_gpu != 0) {
        printf("Application of Q^T to b failed, info: %d\n", info_gpu);
    }

    // solve R*x = Q^T*b
    CHECK_CUBLAS(cublasDtrsm(cublasH, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, 1, &alpha, d_A, m, d_b, n), cublasH);

    // Copy result back 
    CHECK_CUDA(cudaMemcpy(b, d_b, sizeof(double) * n, cudaMemcpyDeviceToHost));

    
    printf("Solution:\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", b[i]);
    }

    // Cleanup
    cudaFree(d_A);
    cudaFree(d_tau);
    cudaFree(d_b);
    cudaFree(devInfo);
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    cublasDestroy(cublasH);

    return 0;
}
