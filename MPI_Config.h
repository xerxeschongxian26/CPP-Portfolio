#ifndef MPI_CONFIG
#define MPI_CONFIG

#include <mpi.h>

/*Declarations for functions found in MPI_Config.cpp;
============================================================================================================================================================================*/ 

bool ProcessCheck(const int &np, const int &Px, const int &Py);
bool cflCheck(double deltat, double Re, double dx, double dy);

void GetNeighbour(MPI_Comm mygrid, int *Neighbour);
void Staffer(const int &rank, const int &Nx, const int &Ny, const int &Px, const int &Py, const double &dx, const double &dy, int *coords, double *start, int &nx, int &ny);


#endif