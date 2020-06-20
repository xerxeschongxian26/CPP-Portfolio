#ifndef LID_DRIVEN_CAVITY
#define LID_DRIVEN_CAVITY

#include <string>
#include "PoissonSolver.h"
#include <mpi.h>

using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity(MPI_Comm mygrid, int rank, int *coords, double *start, int *Neighbour, int nx,
	       	int ny, double deltat, double finalt, double re, double dx, double dy);
    ~LidDrivenCavity();
    
    //Print Functions
    void printVort();
    void printStream();
    
    //Setting Variable Functions
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void SetDiagonals(double dx, double dy);

    //Main and Sub Functions of Solver
    void Initialise();
    void GetBoundaryVec(double *b, char matrix, char x);
    void SetBoundary();
    void SetInteriorVorticity();
    void SetNextInteriorVorticity(); 
    void Integrate();
    void SendRec(double *x, double *x_up, double *x_left,double *x_down, double *x_right);
    void GetVelocity();
    bool GetResidual(double t);
    void FULLSOLVER();


    //Other Functions
    void outputFile(double Lx, double Ly, int Px, int Py);


private:

    MPI_Comm mygrid;                            // New communicator mygrid
    int      rank;
    int      coords[2];       			// 2 dimensional coordinate in Cartesian topology
    double   start[2];     			// Starting coordinates in global coordinate
    int      Neighbour[4];    			// ranks of neighborhood

    double* v       = nullptr;   		// Vorticity function
    double* s       = nullptr;   		// Stream function
    double *s_prev  = nullptr; 			// Store previous s for residual

    double *v_up    = nullptr;
    double *v_down  = nullptr;
    double *v_left  = nullptr;
    double *v_right = nullptr;
    double *s_up    = nullptr;
    double *s_down  = nullptr;
    double *s_left  = nullptr;
    double *s_right = nullptr;


    double *Vx      = nullptr;	 		// Horizontal velocity 
    double *Vy      = nullptr;	 		// Vertical velocity
    
    
    double Lx;  	 			// Length of x-domain
    double Ly;  	 			// Length of y-domain
    int    Nx;  	 			// No. of x-grid points
    int    Ny;  	 			// No. of y-grid points 
    double dt;           			// Time step
    double T;   	 			// Final time
    double Re;  	 			// Reynolds Number

    double dx;		 			// x-direction grid spacing
    double dy;	         			// y-direction grid spacing
    const double U=1.0;   			// Lid velocity
    double tol = 0.001; 		        	// Tolerance criteria of 0.01% 
    double s_nrm= 0.0;				// 2-norm of s
    double residual = 0.0;			// 2-norm of difference between s and s_old 
    double max_residual= 0.0;			// maximum residual along the whole domain

    
    double *A      = nullptr;	
    double *B      = nullptr;    
    double *C      = nullptr;    
    
    double A_MainDiag; 				// A Main Diagonal Constant
    double A_SubDiag;  				// A Sub Diagonal Constant
    double A_SupDiag;   			// A Super Diagonal Constant	
    double B_SubDiag;     			// B Sub Diagonal Constant
    double C_SupDiag;   			// C Super Diagonal Constant  
   
    int    size;				// Size of the grid formed by Nx and Ny	
    int    ldA; 				// Leading dimension of matrix A
    int    ldB; 				// Leading dimension of matrix B
    int    ldC;					// Leading dimension of matrix C

    PoissonSolver *solver_Poisson = nullptr;
};

#endif