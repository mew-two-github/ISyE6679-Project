#ifndef CUSOLVER_STUFF_H

#define CUSOLVER_STUFF_H

// void convertToCSR(double* matrix, int rows, int cols, int* row_ptr, int* col_idx, double* values, int nnz);
extern "C" void solve(double* matrix,  double* b_arr, int m, double* res );

#endif