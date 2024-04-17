#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>


using namespace std;
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
            if (matrix[j * rows + i] != 0) {
                values[k] = matrix[j * rows + i];
                col_idx[k] = j+1;
                k++;
            }
        }
        row_ptr[i + 1] = k+1;
    }
}
int main(int argc, char *argv[]) {
    const int m = 3;
    // const int nnzA = 7;
    int rows = m, cols = 3;
    double matrix[] = {1,4,7,2,5,8,3,6,9};
    int nnzA = 0;
        for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            if (matrix[j * rows + i] != 0) {
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
    
    for(int i=0;i<rows+1;i++)    
    {
        csrRowPtrAcpy[i] = row_ptr[i];
    }
    cout<<nnzA<<"\n";
    for(int i =0;i<nnzA;i++)
    {
        csrValAcpy[i]= values[i];
        csrColIndAcpy[i] = col_idx[i];
    }

    for(int i= 0; i < nnzA; ++i)
    {
        cout<<col_idx[i]<<" "<<values[i]<<"\n";
    }
    cout<<endl;
    for(int i= 0; i < rows+1;++i)
    {
        cout<<row_ptr[i]<<" ";
    }
    cout<<endl;
    return 0;
}
