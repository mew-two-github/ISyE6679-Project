#include "csv.h"
#include "define2.h"
#include <iostream> 
#include <vector>
#include<read_data.h>
#include <codi.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>


using Real = codi::RealForward;
namespace ublas = boost::numeric::ublas;




int main()
{
	
	Real res[1];
	res[0] = 0;
	int Nb = 5, sizeY = 3;
	ublas::compressed_matrix<double> lambda(2*Nb, 1), mu(sizeY, 1), gamma(sizeY, 1), Z(sizeY, 1);
	
	for(size_t i = 0; i < lambda.size1(); ++i) {
        for(size_t j = 0; j < lambda.size2(); ++j) {
             lambda(i,j) = int(i) + int(j);
		}
    }

    Real GX[5];
    for(int i = 0;i<Nb;++i)
    {
    	GX[i] = 29*i;
    }

    size_t j = 0;
    int counter = 0;
	for(size_t i=0 ; i < lambda.size1();i++)
	 {
	 	res[0] += lambda(i,j)*GX[counter];
	 	counter++;
	 }

	 counter = 0;

	Real HX[3] = {225, 154, 3267};
 	for(size_t i = 0; i < mu.size1(); ++i) {
	    for(size_t j = 0; j < mu.size2(); ++j) {
	         mu(i,j) = pow(2,i);
	         gamma(i,j) = pow(3,i);
	         Z(i,j) = 5*i*i - 6*i + 1;
		}
	}
	 for(size_t i=0 ; i < mu.size1();i++)
	 {
	 	res[0] += mu(i,j)*HX[i];
	 	res[0] -= gamma(i,j)*Z(i,j);
	 }
	 cout<<res[0].value()<<endl; // Correct answer is 14389
}