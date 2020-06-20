#ifndef POISSON_SOLVER
#define POISSON_SOLVER

class PoissonSolver
{
public:
	PoissonSolver(const int &Nx, const int &Ny, const double &dx, const double &dy);
	~PoissonSolver();
	
	// Print Functions
	void printVector(int n, double* u);
	void printMatrix(int n, int m, double* M);

	// Cholesky Factorisation Method Functions
	void SetBoundary_Poisson(const double *up, const double *left, const double *down, const double *right);
	void solver_CholFac(double *x, const double *f);
	
private:
	int    Nx;		// No. of x-grid points
	int    Ny;		// No. of y-grid points
	double dx; 		// x-direction grid spacing
	double dy; 		// y-direction grid spacing
	int    size;		// Size of the grid formed by Nx and Ny
	int    ldA;		// Leading dimension of A matrix
	int    info;     

	double *A = nullptr;	// A matrix pointer 
	double *b = nullptr;	// Store boundary vector for cblas routines

	double A_MainDiag;
	double A_SubDiag;	// Sub Diagonal of Matrix A
	double A_SupDiag;	// Super Diagonal of Matrix A

	
};

#endif































