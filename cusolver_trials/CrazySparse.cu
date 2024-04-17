#include <cstdio>
#include <cstdlib>
#include <vector>
#include<iostream>

#include <cuda_runtime.h>
#include <cusolverSp.h>
#include <cusparse.h>

using namespace std;

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

void convertToCSR(double* matrix, int rows, int cols, int* row_ptr, int* col_idx, double* values, int nnz) 
{
    int i, j, k;
   
    /*row_ptr = (int*)malloc((rows + 1) * sizeof(int));
    col_idx = (int*)malloc(nnz * sizeof(int));
    values = (double*)malloc(nnz * sizeof(double));*/

    
    k = 0;
    row_ptr[0] = 1;
    for (i = 0; i < rows; i++) 
    {
        for (j = 0; j < cols; j++) {
            if (abs(matrix[j * rows + i]) <= 1e-8) {
                values[k] = matrix[j * rows + i];
                col_idx[k] = j+1;
                k++;
            }
        }
        row_ptr[i + 1] = k+1;
    }
}

int main(double* matrix,  double* b_arr, int m, ) {
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
    // const int nnzA = 7;
    int rows = m, cols = m;
    
    int nnzA = 0;
        for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            if (abs(matrix[j * rows + i]) <= 1e-8) {
                nnzA++;
            }
        }
    }
    int row_ptr[rows+1];
    int col_idx[nnzA];
    double values[nnzA];


    convertToCSR(matrix, rows, cols, row_ptr, col_idx, values, nnzA);

    std::vector<int> csrRowPtrAcpy(rows+1);
    std::vector<int> csrColIndAcpy(nnzA);
    std::vector<double> csrValAcpy(nnzA);
    // std::vector<double> b = {1,2,3};
    cout<<"Rowptrs\n";
    for(int i=0;i<rows+1;i++)    
    {
        csrRowPtrAcpy[i] = row_ptr[i];
    cout<<csrRowPtrAcpy[i]<<" ";
    }
    cout<<endl;
    cout<<"nnzA"<<nnzA;
    cout<<"\nCol Idx\n";
    for(int i =0;i<nnzA;i++)
    {
        csrValAcpy[i]= values[i];
        csrColIndAcpy[i] = col_idx[i];
    cout<<csrColIndAcpy[i]<<" ";
    }
    cout<<endl;
    const std::vector<int> csrRowPtrA = csrRowPtrAcpy;
    const std::vector<int> csrColIndA = csrColIndAcpy;
    const std::vector<double> csrValA = csrValAcpy;
    // const std::vector<double> b = {40,40,40,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,225.753,13.4936,-96.244,-54.9104,-7.6,5.38456,-2.84615e-311,24,-35.3885,-11.0896,-3.5,-8.85719,-14.7368,-14.9,25.3857,14.4312,-8.28304,29.5385,-1.66932,-1.52778,10.9658,-1.8,13.2174,-0.908676,-5,2.49465,2.51826,-5.01828e-310,0,0,0,0,0,0,0,0,0,0,0,0};
    vector<double> b(m);
    for(int i = 0;i < m;++i)
    {
        b[i] = b_arr[i];
    }
    const int batchSize = 17;
    cout<<"line 126-130"<<endl;


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
    cout<<"line 135"<<endl;
    batchSizeMax = std::min(batchSizeMax, batchSize);
    batchSizeMax = 2;

    CUSOLVER_CHECK(cusolverSpDcsrqrBufferInfoBatched(cusolverH, m, m, nnzA, descrA, d_csrValA, d_csrRowPtrA,
                                                     d_csrColIndA, batchSizeMax, info, &size_internal, &size_qr));

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&buffer_qr), size_qr));
    cout<<"line 142"<<endl;
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

int main()
{
    double matrix[] = {0.0860585,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,1,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,1,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6.27048,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,20.2695,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5.29894,-4.5994,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,16.1788,15.3409,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,-1.18607,-2.02373,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,4.99719,5.16531,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,-1.76187,-2.00586,-5.04306,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,5.34622,5.11969,28.9878,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.08756,-1.77754,0,-5.03622,-1.43187e-310,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,4.4891,5.42773,0,21.989,4.24576,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,-5.23268e-312,-6.93031,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,4.04736,14.0862,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,4.87286,0,-6.65675e-310,15.787,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,6.02897,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-20.3937,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,1.83216,0,-6.65675e-310,9.65358,0,8.14312,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-18.9379,-1.98822,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,5.0569,4.65397,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-2.09185,10.1942,0,-14.8174,-1.97693,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,4.38058,1.73056e-311,0,-5.88854,4.62755,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.63282,10.1942,0,-14.8174,3.4011e-309,0,-2.61345,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,3.39832,1.73056e-311,0,-5.88854,1.45298e-309,0,2.3646,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-3.31593,10.1942,0,-14.8174,3.4011e-309,0,-2.62589,-1.17793,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,6.53021,1.73056e-311,0,-5.88854,1.45298e-309,0,2.37586,2.39834,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-16.3212,3.4011e-309,0,-1.58101e-321,-1.19385,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-2.68991,1.45298e-309,0,6.37345e-322,2.43075,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,6.37345e-322,2.39223e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.30357e-312,3.99431e-312,0,1.93476,2.80903e-309,-1.32668e-309,10.1942,0,-14.8174,3.4011e-309,0,-1.58101e-321,4.87072e-310,0,8.00344e-313,1.30811e-312,0,5.13e-312,0,-6.65675e-310,1.73056e-311,0,-5.88854,1.45298e-309,0,8.63281e-320,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,-6.27048,-5.29894,3.30357e-312,3.30357e-312,-1.08756,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,3.30357e-312,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.99431e-312,1,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,-4.5994,-1.18607,-1.76187,-1.77754,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,3.99431e-312,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-2.02373,-2.00586,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,-5.04306,-5.03622,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,1.93476,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,-1.43187e-310,-5.23268e-312,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,2.80903e-309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-6.93031,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-2.09185,-1.63282,-3.31593,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,-1.32668e-309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,10.1942,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-20.3937,-18.9379,-14.8174,-14.8174,-14.8174,-16.3212,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,-14.8174,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,-1.98822,-1.97693,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,3.4011e-309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-2.61345,-2.62589,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,-1.58101e-321,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,-1.17793,-1.19385,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,4.87072e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,20.2695,16.1788,8.00344e-313,8.00344e-313,4.4891,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,8.00344e-313,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,15.3409,4.99719,5.34622,5.42773,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,1.30811e-312,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.16531,5.11969,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,28.9878,21.989,5.13e-312,4.87286,5.13e-312,1.83216,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,5.13e-312,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.24576,4.04736,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,14.0862,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,4.38058,3.39832,6.53021,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,-6.65675e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,15.787,6.02897,9.65358,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,1.73056e-311,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,8.14312,5.0569,-5.88854,-5.88854,-5.88854,-2.68991,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,-5.88854,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,4.65397,4.62755,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,1.45298e-309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8.63281e-320,6.37345e-322,6.37345e-322,8.63281e-320,6.37345e-322,6.37345e-322,8.63281e-320,6.37345e-322,6.37345e-322,8.63281e-320,6.37345e-322,6.37345e-322,8.63281e-320,6.37345e-322,6.37345e-322,8.63281e-320,6.37345e-322,2.3646,2.37586,6.37345e-322,6.37345e-322,8.63281e-320,6.37345e-322,6.37345e-322,8.63281e-320,6.37345e-322,6.37345e-322,8.63281e-320,6.37345e-322,6.37345e-322,8.63281e-320,6.37345e-322,6.37345e-322,8.63281e-320,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8.88581e-310,2.39223e-310,2.39223e-310,8.88581e-310,2.39223e-310,2.39223e-310,8.88581e-310,2.39223e-310,2.39223e-310,8.88581e-310,2.39223e-310,2.39223e-310,8.88581e-310,2.39223e-310,2.39223e-310,8.88581e-310,2.39223e-310,2.39223e-310,2.39834,2.43075,2.39223e-310,8.88581e-310,2.39223e-310,2.39223e-310,8.88581e-310,2.39223e-310,2.39223e-310,8.88581e-310,2.39223e-310,2.39223e-310,8.88581e-310,2.39223e-310,2.39223e-310,8.88581e-310,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double b[] = {40,40,40,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,225.753,13.4936,-96.244,-54.9104,-7.6,5.38456,-2.84615e-311,24,-35.3885,-11.0896,-3.5,-8.85719,-14.7368,-14.9,25.3857,14.4312,-8.28304,29.5385,-1.66932,-1.52778,10.9658,-1.8,13.2174,-0.908676,-5,2.49465,2.51826,-5.01828e-310,0,0,0,0,0,0,0,0,0,0,0,0}; 
}