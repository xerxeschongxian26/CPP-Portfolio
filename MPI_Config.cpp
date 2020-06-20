#include <iostream>
#include <mpi.h>
#include "MPI_Config.h"

using namespace std;

/*Functions to perform process and CFL check;
================================================================================================================================ */

// No. of processors is compatible with no. of partitions.
bool ProcessCheck(const int &np, const int &Px, const int &Py)
{
	if (np == Px * Py) {
        return true;
    }
	else{
        return false;
    }
}

// Input timestep satisfies CFL condition
bool cflCheck(double deltat, double Re, double dx, double dy) 		
{
    if (deltat<(Re/4)*dx*dy) {
        return true;
    }
    else {
        return false;
    }
}

/*Use of Cart_shift to obtain surrounding 4 ranks of each process/rank, stored in elements of Neighbour;
================================================================================================================================ */

void GetNeighbour(MPI_Comm mygrid, int *Neighbour)
{
	MPI_Cart_shift(mygrid, 0, -1, Neighbour ,    Neighbour + 2);				//Direction 0: Row, top[0] bottom[2]
	MPI_Cart_shift(mygrid, 1,  1, Neighbour + 1, Neighbour + 3);				//Direction 1: Col, left[1] right[3]
}

/*Function to distribute work to each rank.;
================================================================================================================================ */
void Staffer(const int &rank, const int &Nx, const int &Ny, const int &Px, const int &Py,
            const double &dx, const double &dy, int *coords, double *start, int &nx, int &ny)
{
    // ranks with smaller coordinates have 1 larger size in grids
    nx = (Nx-2) / Px;
    ny = (Ny-2) / Py;
    int remainder_x = (Nx-2) % Px;
    int remainder_y = (Ny-2) % Py;

    // Caution: start[i, j] and coords[j, i]. This is due to the axis defined in the problem.  
    // x direction 
    if (coords[1] < remainder_x){
        nx++;
        start[0] = (coords[1]*nx + 1)*dx;
    }
    else{
        start[0] = (Nx-1)*dx - (Px - coords[1])*dx*nx;
    }
    // y direction
    if (coords[0] < remainder_y){
        ny++;
        start[1] = (coords[0]*ny + 1)*dy;
    }
    else
    {
        start[1] = (Ny-1)*dy - (Py - coords[0])*dy*ny;
    } 
}



