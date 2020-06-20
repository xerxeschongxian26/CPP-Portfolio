#include<cstring>
#include<iostream>
#include<iomanip>
#include<cmath>

#include "PoissonSolver.h"
#include "cblas.h"

using namespace std;
 
#define F77NAME(x) x##_
extern "C"
{
	void F77NAME(dpbtrf) (const char &UPLO, const int &N, const int &KD, const double *AB,
		       const int &LDAB,	int &INFO);
	void F77NAME(dpbtrs) (const char &UPLO, const int &N, const int &KD, const int &NRHS,
		       	const double *AB, const int &LDAB, double *B, const int &LDB, int &INFO);
}

PoissonSolver::PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy)
{
	this->Nx = Nx;
	this->Ny = Ny;
	this->dx = dx;
	this->dy = dy;

	// Get Matrix A
	ldA = Ny + 1;
	size = Nx * Ny;
	A = new double[size * ldA];
	cblas_dscal(size*ldA, 0.0, A, 1);	
	A_MainDiag = 2/(dx*dx) + 2/(dy*dy);
	A_SubDiag = -1/(dy*dy);	
	A_SupDiag = -1/(dx*dx);	
	for (int j = 0; j < size; j++)
	{
		A[ldA-1 + j*ldA] = A_MainDiag;
		if (j % Ny != 0)   {
			A[ldA-2 + j*ldA] = A_SubDiag;
		}
		if (j >= Ny){
			A[j * ldA] = A_SupDiag;
		}
	}
	info = 1;
}
PoissonSolver::~PoissonSolver()
{
	delete[] A;
	delete[] b;
}

/*Print Functions;
================================================================================================================================ */

void PoissonSolver::printVector(int n, double* u) {
   for (int i = 0; i < n; ++i) {
       cout << u[i] << endl;
   }
   cout << endl;
}
void PoissonSolver::printMatrix(int n, int m, double* M) {
   for (int j = 0; j < m; ++j) {
       for (int i = 0; i < n; ++i) {
           cout << setw(6) << M[j+i*m] << " ";
       }
       cout << endl;
   }
}

/*Cholesky Factorisation Method Functions;
================================================================================================================================ */

// Creates a boundary vector b from conditions of 4 neighbours (up,left,down,right)
void PoissonSolver::SetBoundary_Poisson(const double *up, const double *left, const double *down, const double *right)
{	
	b = new double[size];
	cblas_dscal(size, 0.0, b, 1);
	for (int i = 0; i < Nx; i++)
	{
		b[i*Ny] += down[i] * A_SubDiag;			
		b[i*Ny + Ny-1] += up[i] * A_SubDiag; 			
	}
	for (int i = 0; i < Ny; i++)
	{
		b[i] += left[i] * A_SupDiag; 				
		b[i+(Nx-1)*Ny] += right[i] * A_SupDiag;  	
	}
}

//Solve using the Cholesky Factorisation Method
void PoissonSolver::solver_CholFac(double *x, const double *f)
{
	// Check if factorisation occured only once
	if (info != 0)	
	{
		F77NAME(dpbtrf) ('U', size , Ny, A, ldA, info);
	}
	cblas_daxpy(size, -1.0, f, 1, b, 1);
	cblas_dscal(size, -1.0, b, 1);

	F77NAME(dpbtrs) ('U', size, Ny, 1, A, ldA, b, size, info);
	cblas_dcopy(size, b, 1, x, 1); 							//Output x is used to update s
	
}
