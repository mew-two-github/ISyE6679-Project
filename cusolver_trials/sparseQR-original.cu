#include <cstdio>
#include <cstdlib>
#include <vector>

#include <cuda_runtime.h>
#include <cusolverSp.h>
#include <cusparse.h>

// Error handling macros
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            printf("CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

#include <iostream>
#define CUSOLVER_CHECK(call) { \
    cusolverStatus_t err; \
    if ((err = (call)) != CUSOLVER_STATUS_SUCCESS) { \
        std::cerr << "CUSOLVER error: " << err << ", " << "file: " << __FILE__ << ", line: " << __LINE__ << std::endl; \
        exit(1); \
    } \
}

#define CUSPARSE_CHECK(call) \
    do { \
        cusparseStatus_t error = call; \
        if (error != CUSPARSE_STATUS_SUCCESS) { \
            printf("CUSPARSE error at %s:%d\n", __FILE__, __LINE__); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

int main(int argc, char *argv[]) {
    cusolverSpHandle_t cusolverH = NULL;
    csrqrInfo_t info = NULL;
    cusparseMatDescr_t descrA = NULL;
    cudaStream_t stream = NULL;

    int *d_csrRowPtrA = nullptr;
    int *d_csrColIndA = nullptr;
    double *d_csrValA = nullptr;
    double *d_b = nullptr; // batchSize * m
    double *d_x = nullptr; // batchSize * m

    size_t size_qr = 0;
    size_t size_internal = 0;
    void *buffer_qr = nullptr; // working space for numerical factorization


    const int m = 4;
    const int nnzA = 7;
    const std::vector<int> csrRowPtrA = {1, 2, 3, 4, 8};
    const std::vector<int> csrColIndA = {1, 2, 3, 1, 2, 3, 4};
    const std::vector<double> csrValA = {1.0, 2.0, 3.0, 0.1, 0.1, 0.1, 4.0};
    const std::vector<double> b = {1.0, 1.0, 1.0, 1.0};
    
    
    const int batchSize = 17;



    std::vector<double> csrValABatch(nnzA * batchSize);
    std::vector<double> bBatch(m * batchSize);
    std::vector<double> xBatch(m * batchSize);

    for (int colidx = 0; colidx < nnzA; colidx++) {
        double Areg = csrValA[colidx];
        for (int batchId = 0; batchId < batchSize; batchId++) {
            double eps = (static_cast<double>((std::rand() % 100) + 1)) * 1.e-4;
            csrValABatch[batchId * nnzA + colidx] = Areg + eps;
        }
    }

    for (int j = 0; j < m; j++) {
        double breg = b[j];
        for (int batchId = 0; batchId < batchSize; batchId++) {
            double eps = (static_cast<double>((std::rand() % 100) + 1)) * 1.e-4;
            bBatch[batchId * m + j] = breg + eps;
        }
    }

    CUSOLVER_CHECK(cusolverSpCreate(&cusolverH));
    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    CUSOLVER_CHECK(cusolverSpSetStream(cusolverH, stream));
    CUSPARSE_CHECK(cusparseCreateMatDescr(&descrA));
    CUSPARSE_CHECK(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
    CUSPARSE_CHECK(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ONE));
    CUSOLVER_CHECK(cusolverSpCreateCsrqrInfo(&info));

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_csrValA), sizeof(double) * csrValABatch.size()));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_csrColIndA), sizeof(int) * csrColIndA.size()));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_csrRowPtrA), sizeof(int) * csrRowPtrA.size()));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_b), sizeof(double) * bBatch.size()));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_x), sizeof(double) * xBatch.size()));

    CUDA_CHECK(cudaMemcpyAsync(d_csrColIndA, csrColIndA.data(), sizeof(int) * csrColIndA.size(),
                               cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_csrRowPtrA, csrRowPtrA.data(), sizeof(int) * csrRowPtrA.size(),
                               cudaMemcpyHostToDevice, stream));

    CUSOLVER_CHECK(cusolverSpXcsrqrAnalysisBatched(cusolverH, m, m, nnzA, descrA, d_csrRowPtrA,
                                                   d_csrColIndA, info));

    CUDA_CHECK(cudaStreamSynchronize(stream));

    size_t free_mem = 0;
    size_t total_mem = 0;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));

    int batchSizeMax = 2;
    while (batchSizeMax < batchSize) {
        CUSOLVER_CHECK(cusolverSpDcsrqrBufferInfoBatched(cusolverH, m, m, nnzA, descrA, d_csrValA, d_csrRowPtrA,
                                                         d_csrColIndA, batchSizeMax, info, &size_internal, &size_qr));
        if ((size_internal + size_qr) > free_mem) {
            batchSizeMax /= 2;
            break;
        }
        batchSizeMax *= 2;
    }

    batchSizeMax = std::min(batchSizeMax, batchSize);
    batchSizeMax = 2;

    CUSOLVER_CHECK(cusolverSpDcsrqrBufferInfoBatched(cusolverH, m, m, nnzA, descrA, d_csrValA, d_csrRowPtrA,
                                                     d_csrColIndA, batchSizeMax, info, &size_internal, &size_qr));

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&buffer_qr), size_qr));

    for (int idx = 0; idx < batchSize; idx += batchSizeMax) {
        const int cur_batchSize = std::min(batchSizeMax, batchSize - idx);
        CUDA_CHECK(cudaMemcpyAsync(d_csrValA, csrValABatch.data() + idx * nnzA,
                                   sizeof(double) * nnzA * cur_batchSize, cudaMemcpyHostToDevice,
                                   stream));
        CUDA_CHECK(cudaMemcpyAsync(d_b, bBatch.data() + idx * m, sizeof(double) * m * cur_batchSize,
                                   cudaMemcpyHostToDevice, stream));

        CUSOLVER_CHECK(cusolverSpDcsrqrsvBatched(cusolverH, m, m, nnzA, descrA, d_csrValA,
                                                 d_csrRowPtrA, d_csrColIndA, d_b, d_x,
                                                 cur_batchSize, info, buffer_qr));

        CUDA_CHECK(cudaMemcpyAsync(xBatch.data() + idx * m, d_x, sizeof(double) * m * cur_batchSize,
                                   cudaMemcpyDeviceToHost, stream));
    }

    CUDA_CHECK(cudaStreamSynchronize(stream));

    const int baseA = (CUSPARSE_INDEX_BASE_ONE == cusparseGetMatIndexBase(descrA)) ? 1 : 0;

    for (int batchId = 0; batchId < batchSize; batchId++) {
        double *csrValAj = csrValABatch.data() + batchId * nnzA;
        double *xj = xBatch.data() + batchId * m;
        double *bj = bBatch.data() + batchId * m;
        double sup_res = 0;
        for (int row = 0; row < m; row++) {
            const int start = csrRowPtrA[row] - baseA;
            const int end = csrRowPtrA[row + 1] - baseA;
            double Ax = 0.0;
            for (int colidx = start; colidx < end; colidx++) {
                const int col = csrColIndA[colidx] - baseA;
                const double Areg = csrValAj[colidx];
                const double xreg = xj[col];
                Ax += Areg * xreg;
            }
            double r = bj[row] - Ax;
            sup_res = std::max(sup_res, fabs(r));
        }
        printf("batchId %d: sup|bj - Aj*xj| = %E \n", batchId, sup_res);
    }

    for (int batchId = 0; batchId < batchSize; batchId++) {
        double *xj = xBatch.data() + batchId * m;
        for (int row = 0; row < m; row++) {
            printf("x%d[%d] = %E\n", batchId, row, xj[row]);
        }
        printf("\n");
    }

    CUDA_CHECK(cudaFree(d_csrRowPtrA));
    CUDA_CHECK(cudaFree(d_csrColIndA));
    CUDA_CHECK(cudaFree(d_csrValA));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(buffer_qr));

    CUSOLVER_CHECK(cusolverSpDestroy(cusolverH));
    CUDA_CHECK(cudaStreamDestroy(stream));
    CUDA_CHECK(cudaDeviceReset());

    return EXIT_SUCCESS;
}
