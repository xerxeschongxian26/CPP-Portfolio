#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include "cblas.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <mpi.h>

using namespace std;

#define FORMAT(vOut, x, y, v, s, Vx, Vy) vOut << setw(10) << x << setw(10) << y << setw(12) << v << setw(12) << s << setw(12) << Vx << setw(12) << Vy << endl;

LidDrivenCavity::LidDrivenCavity(MPI_Comm mygrid, int rank, int *coords, double *start, int *Neighbour, int nx,	int ny, double deltat, double finalt, double re, double dx, double dy)
{   this-> mygrid = mygrid;
    this-> rank = rank;
    memcpy(this->coords, coords, 2*sizeof(int));
    memcpy(this->Neighbour, Neighbour, 4*sizeof(int));
    memcpy(this->start, start, 2*sizeof(double));
    Nx = nx;
    Ny = ny;
    size = Nx * Ny;
    ldA = 1+Ny;
    ldB = 3;
    ldC = 1+2*Ny;
    T = finalt;
    Re = re;
    this->dx = dx;
    this->dy = dy;

    SetDomainSize(Lx, Ly);
    SetTimeStep(deltat);
    SetFinalTime(T);
    SetReynoldsNumber(Re);
    SetDiagonals(dx, dy);
}
LidDrivenCavity::~LidDrivenCavity()
{
	delete[] v;
	delete[] s;
	delete[] v_up;
	delete[] v_down;
	delete[] v_left;
	delete[] v_right;
	delete[] s_up;
	delete[] s_down;
	delete[] s_left;
	delete[] s_right;
	delete[] A;
	delete[] B;
	delete[] C;
}

/*Print Functions;
================================================================================================================================ */

//Used to print matrix, for checking and displaying solution.
void LidDrivenCavity::printVort() {
    cout.precision(4);
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            cout << setw(6) << v[j*Nx+i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
//Used to print matrix, for checking and displaying solution.
void LidDrivenCavity::printStream() {
    cout.precision(4);
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j <Ny; ++j) {
            cout << setw(6) << s[j*Nx+i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

/*Setting Variable Functions;
================================================================================================================================ */

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx = xlen;
    Ly = ylen;
}
void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    Nx=nx;
    Ny=ny;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    dt = deltat;
}
void LidDrivenCavity::SetFinalTime(double finalt)
{
    T= finalt;
}
void LidDrivenCavity::SetReynoldsNumber(double re)
{
    Re=re;
}
void LidDrivenCavity::SetDiagonals(double dx, double dy)
{
   A_MainDiag  = 2/(dx*dx) + 2/(dy*dy);			
   A_SubDiag  = -1/(dy*dy);					
   A_SupDiag   = -1/(dx*dx);						
   B_SubDiag  = 1/(2*dy);
   C_SupDiag   = 1/(2*dx);		
}

/*Main and Sub Functions of Solver;
================================================================================================================================ */

void LidDrivenCavity::Initialise()
{
// Allocate the dynamics memory for vorticity,stream and the matrices need to solve the system. Scale by zero
    v      = new double[size];
    s	   = new double[size];
    s_prev = new double[size];

    Vx    = new double[size];						//Velocities in x-direction
    Vy    = new double[size];						//Velocities in y-direction

    A     = new double[size * ldA];					// Store as symmetric tri-diagonal banded matrix with bwidth = Ny
    B     = new double[size * ldB];    					// Store as banded matrix with bwidth = 1
    C     = new double[size * ldC];					// Store as banded matrix with bwidth = Ny

    cblas_dscal(size, 0.0, v, 1);					//Scale by zero
    cblas_dscal(size, 0.0, s, 1);					//Scale by zero 

// Allocate the dynamic memory for BCs  of v and s. Scale by zero
    v_up     = new double[Nx];
    v_down   = new double[Nx];
    v_left   = new double[Ny];
    v_right  = new double[Ny];
    s_up     = new double[Nx];
    s_down   = new double[Nx];
    s_left   = new double[Ny];
    s_right  = new double[Ny];

    cblas_dscal(Nx, 0.0, v_up, 1);
    cblas_dscal(Nx, 0.0, v_down, 1);
    cblas_dscal(Ny, 0.0, v_left, 1);
    cblas_dscal(Ny, 0.0, v_right, 1);
    cblas_dscal(Nx, 0.0, s_up, 1);
    cblas_dscal(Nx, 0.0, s_down, 1);
    cblas_dscal(Ny, 0.0, s_left, 1);
    cblas_dscal(Ny, 0.0, s_right, 1);

// Populate matrices A B C with their known diagonal constants
    for (int i = 0; i < size; i++){
	A[ldA-1 + i*ldA] = A_MainDiag;
	if (i % Ny != 0){
		A[ldA-2 + i*ldA] = A_SubDiag;
	}
	if (i >= Ny){
		A[i * ldA] = A_SupDiag;
	}
    }
    for (int i = 0; i < size; i++){
	if (i % Ny != 0){
	    B[i*ldB] = B_SubDiag;
	}
	if (i % Ny != Ny-1){
	    B[2 + i*ldB] = -B_SubDiag;
	}
    }
    for (int i = 0; i < size; i++){
	if (i >= Ny){
	    C[i*ldC] = C_SupDiag;
	}
	if (i < size-Ny){   
	    C[ldC-1 + i*ldC] = -C_SupDiag;
	}
    }
}

//Function to output boundary vector b in the linear system y = Ax + b. Input matrix can be either A, B or C.
// Char value ('v' 's') indicates if is the vorticity or stream function
void LidDrivenCavity::GetBoundaryVec(double *b, char matrix, char x) 
{
	switch(matrix)
	{
		case 'A':
		{
			if (x == 's')					//Streamfunction
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] += s_left[i] * A_SupDiag;
					b[i+(Nx-1)*Ny] += s_right[i] * A_SupDiag;
				}
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] += s_down[i] * A_SubDiag;
					b[i*Ny + Ny-1] += s_up[i] * A_SubDiag;
				}
			}
			else if (x == 'v')				//Vorticity
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] += v_left[i] * A_SupDiag;
					b[i+(Nx-1)*Ny] += v_right[i] * A_SupDiag;
				}
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] += v_down[i] * A_SubDiag;
					b[i*Ny + Ny-1] += v_up[i] * A_SubDiag;
				}
			} break;
		}
		case 'B':
		{
			if (x == 's')
			{
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] = s_down[i] * (-B_SubDiag);
					b[i*Ny + Ny-1] = s_up[i] * B_SubDiag;
				}
			}
			else if (x == 'v')
			{
				for (int i = 0; i < Nx; i++)
				{
					b[i*Ny] = v_down[i] * (-B_SubDiag);
					b[i*Ny + Ny-1] = v_up[i] * B_SubDiag;
				}
			} break;
		}
		case 'C':
		{
			if (x == 's')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] = s_left[i] * (-C_SupDiag);
					b[i+(Nx-1)*Ny] = s_right[i] * C_SupDiag;
				}
			}
			else if (x == 'v')
			{
				for (int i = 0; i < Ny; i++)
				{
					b[i] = v_left[i] * (-C_SupDiag);
					b[i+(Nx-1)*Ny] = v_right[i] * C_SupDiag;
				}
			} break;
		}
		default: break;
	}
}

//Populate v matrix with BCs by screening for ranks with out-of-boundary neighbours (a negative rank)
void LidDrivenCavity:: SetBoundary(){

    if (Neighbour[0] == -2){
        for (int i = 0; i < Nx; i++){
            v_up[i]=(s_up[i]-s[i*Ny+(Ny-1)])*2/(dy*dy)-2*U/dy;	
        }
    }

    if (Neighbour[2] == -2){
        for (int i = 0; i < Nx; i++){
            v_down[i]=(s_down[i]-s[i*Ny])*2/(dy*dy);			
        }
    }

    if (Neighbour[1] == -2){
        for (int i = 0; i < Ny; i++){
            v_left[i] =(s_left[i]-s[i])*2/(dx*dx); 			
        }
    }

    if (Neighbour[3] == -2){
        for (int i = 0; i < Ny; i++){
            v_right[i] =(s_right[i]-s[(Nx-1)*Ny+i])*2/(dx*dx);		
        }
    }
}
//Solves for the interior v values at current time step
void LidDrivenCavity:: SetInteriorVorticity(){
    	double *v_temp =  new double[size];
	GetBoundaryVec(v_temp, 'A', 's');						
	cblas_dsbmv (CblasColMajor, CblasUpper, size, Ny, 1.0, A, ldA, s, 1, 1.0, v_temp, 1); 	 //v_temp=As+v_temp
	cblas_dcopy (size, v_temp, 1, v, 1);							 //v = v_temp
	delete[] v_temp;
}

// Solves for the interior v values at next time step
// 3 distinct products in equation to obtain the interior vorticity at the next time step
// Each product (b_1,b_2,b_3) is evaluated individually and combined using the cblas_daxpy routine
void LidDrivenCavity:: SetNextInteriorVorticity() 
{
	double *b_1 = new double[size];	
	double *b_2 = new double[size];	
	double *b_3 = new double[size];	

	cblas_dscal(size, 0.0, b_1, 1);
	cblas_dscal(size, 0.0, b_2, 1);
	cblas_dscal(size, 0.0, b_3, 1);

	// Calculate b_1
	GetBoundaryVec(b_1, 'A', 'v');
	cblas_dsbmv (CblasColMajor, CblasUpper, size, Ny, 1.0, A, ldA, v, 1, 1.0, b_1, 1);

	// Calculate b_2
	double *temp_1 = new double[size];
	double *temp_2 = new double[size];	
	cblas_dscal(size, 0.0,temp_1,1);
	cblas_dscal(size, 0.0,temp_2,1);

	GetBoundaryVec(temp_1, 'C', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, Ny, Ny, 1.0, C, ldC, s, 1, 1.0, temp_1, 1);

	GetBoundaryVec(temp_2, 'B', 'v');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, 1, 1, 1.0, B, ldB, v, 1, 1.0, temp_2, 1);

	for (int i = 0; i < size; i++)
	{
		b_2[i] = temp_1[i] * temp_2[i];
	}

	// Calculate b_3
	double *temp_3 = new double[size];	
	double *temp_4 = new double[size];	
	cblas_dscal(size,0.0,temp_3,1);
	cblas_dscal(size,0.0,temp_4,1);

	GetBoundaryVec(temp_3, 'C', 'v');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, Ny, Ny, 1.0, C, ldC, v, 1, 1.0, temp_3, 1);

	GetBoundaryVec(temp_4, 'B', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, 1, 1, 1.0, B, ldB, s, 1, 1.0, temp_4, 1);

	for (int i = 0; i < size; i++)
	{
		b_3[i] = temp_3[i] * temp_4[i];
	}
	
	// Calculate interior v values at next time step
	cblas_daxpy (size, -dt/Re, b_1, 1, v, 1);			//Output stored in v
	cblas_daxpy (size,  dt,    b_2, 1, v, 1);			//Output stored in v
	cblas_daxpy (size, -dt,    b_3, 1, v, 1);			//Final output after 3 cblas routines is  v_next = v_current -(dt/Re)*b_1 + dt*b_2 - dt*b_3
	
	//Clear memory
	delete[] b_1;
	delete[] b_2;
	delete[] b_3;
	delete[] temp_1;
	delete[] temp_2;
	delete[] temp_3;
	delete[] temp_4;

}
//Solves poisson problem
void LidDrivenCavity::Integrate(){
	if (solver_Poisson == nullptr)  
	{
		solver_Poisson = new PoissonSolver(Nx, Ny, dx, dy);
	}
	solver_Poisson->SetBoundary_Poisson(s_up, s_left, s_down, s_right);
	solver_Poisson->solver_CholFac(s, v);
}

//Control the flow of information (input/output) in the vertical and horizontal direction of each rank
void LidDrivenCavity::SendRec(double *x, double *x_up, double *x_left, double *x_down, double *x_right)
{
//Vertical Direction
    // U -> D: Source= Neighbour[0], Destination= Neighbour[2]
    double temp_down[Nx];
    cblas_dcopy(Nx, x, Ny, temp_down, 1);
    MPI_Sendrecv(temp_down, Nx, MPI_DOUBLE, Neighbour[2], 2, x_up, Nx, MPI_DOUBLE, Neighbour[0], 2, mygrid, MPI_STATUS_IGNORE);

    //  D -> U: Source= Neighbour[2], Destination= Neighbour[0]
    double temp_up[Nx];
    cblas_dcopy(Nx, x + Ny-1, Ny, temp_up, 1);
    MPI_Sendrecv(temp_up, Nx, MPI_DOUBLE, Neighbour[0], 0, x_down, Nx, MPI_DOUBLE, Neighbour[2], 0, mygrid, MPI_STATUS_IGNORE);

//Horizontal Direction
    // L -> R: Source= Neighbour[1], Destination= Neighbour[3]
    MPI_Sendrecv(x + (Nx-1)*Ny, Ny, MPI_DOUBLE, Neighbour[3], 3, x_left, Ny, MPI_DOUBLE, Neighbour[1], 3, mygrid, MPI_STATUS_IGNORE);

    //  R -> L: Source= Neighbour[3], Destination= Neighbour[1]
    MPI_Sendrecv(x, Ny, MPI_DOUBLE, Neighbour[1], 1, x_right, Ny, MPI_DOUBLE, Neighbour[3], 1, mygrid, MPI_STATUS_IGNORE);


}
 //Get velocities from s
void LidDrivenCavity::GetVelocity()
{
//Vertical Direction
	GetBoundaryVec(Vy, 'C', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, Ny, Ny, 1.0, C, ldC, s, 1, 1.0, Vy, 1);

//Horizontal Direction
	GetBoundaryVec(Vx, 'B', 's');
	cblas_dgbmv (CblasColMajor, CblasNoTrans, size, size, 1, 1, 1.0, B, ldB, s, 1, 1.0, Vx, 1);
}

bool LidDrivenCavity::GetResidual(double t)
{
	s_nrm = cblas_dnrm2(size, s, 1);
        cblas_daxpy(size, -1.0, s, 1, s_prev, 1);						// Difference between prev and current s
        residual = cblas_dnrm2(size, s_prev, 1); 						// 2-norm using cblas routine
        residual /= s_nrm;									// Definition of residual

	MPI_Reduce(&residual, &max_residual, 1, MPI_DOUBLE, MPI_MAX, 0, mygrid);		//Reduces residual variable from all ranks and returns the max value to max_residual located in process 0
        if (rank == 0)
        {
                cout << "At time = " << setw(6) << t << setw(30) << "Max Residual = " << max_residual << endl;
        }
        MPI_Bcast(&max_residual, 1, MPI_DOUBLE, 0, mygrid);     				//Process 0 broadcasts max_residual to all other processes, other process receive the broadcast

        if(max_residual < tol)
        {
	    return 0;                                                            
	}
	else {
	    return 1;
	}
}

void LidDrivenCavity::FULLSOLVER(){

    double t=0.0;										//Initialise time
    bool status_residual;									//Initialise status residual
    Initialise();										//Initialise v,s,s_prev,A,B,C,

    while (t<T)
    {
        SetBoundary();										//Apply BCs
        SetInteriorVorticity();									//Set interior v at current t
        SendRec(v, v_up, v_left, v_down, v_right);		
        SetNextInteriorVorticity();								//Set interior v at next t
	
	cblas_dcopy(size, s, 1, s_prev, 1);							//store streamfunction of the previous step in s_prev
        for (unsigned int i=0;i<5;i++)   	
	{
            Integrate();
            SendRec(s, s_up, s_left, s_down, s_right);
        }
	status_residual=GetResidual(t);
	if (status_residual==0)
	{	 
		if (rank == 0)									//Message output only once at process 0
		{										
			cout << "Tolerance criteria reached: Steady State Solution has been achieved" << endl;
	    	}
		break;					
       	}
        t+=dt;
    }
    GetVelocity();
}

/*Other Functions;
================================================================================================================================ */
// Outputs steady-state solution to a text file
// Preprocessor directive FORMAT not crucial to output of data. It defines a format of outputting data to the file to ensure neatness of outputFile function.
void LidDrivenCavity::outputFile(double Lx, double Ly, int Px, int Py)
{
	for (int k = 0; k < Px *Py; k++)										// Moves through each partition of the domain
	{
		if (k == rank)												// At each rank, begin the data output specified below
		{
			if (k == 0)
			{
				ofstream vOut("data.txt", ios::out | ios::trunc);					// File is open for output ONLY & 'refreshes' the file by deleting before reading data
				FORMAT(vOut, "x", "y", "v", "s", "Vx", "Vy");						// Column headers
				FORMAT(vOut, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);						// Bottom-left corner
				FORMAT(vOut, 0.0, Ly, 0.0, 0.0, 0.0, 0.0);						// Top-left corner
				FORMAT(vOut, Lx, 0.0, 0.0, 0.0, 0.0, 0.0);						// Bottom-right Corner
				FORMAT(vOut, Lx, Ly, 0.0, 0.0, 0.0, 0.0);						// Top-right Corner
				vOut.close();
			}
			ofstream vOut("data.txt", ios::out | ios::app);							// Defines files access modes to "Output.txt"
			double x, y;											// File is open for output ONLY & data is appended to end of file
			
			// Interior
			vOut.precision(5);
			for (int i = 0; i < Nx; i++)
			{
				for (int j = 0; j < Ny; j++)
				{
					x = start[0] + i*dx;								// Array of x-coordinates
					y = start[1] + j*dy;								// Array of y-coordinates
					FORMAT(vOut, x, y, v[i*Ny + j], s[i*Ny + j], Vx[i*Ny + j], Vy[i*Ny + j]);
				}
			}

			// Top: Located at the 'top-most' edge. Correspond to a coords[] row index = number of partitions in y-direction
			if (coords[0] == Py - 1)
			{
				y = Ly;
				for (int i = 0; i < Nx; i++)
				{
					x = start[0] + i*dx;
					FORMAT(vOut, x, y, v_up[i], 0.0, U, 0.0);
				}
			}
			// Bottom: Located at the 'bottom-most' edge. Correspond to a coords[] row index = 0
			if (coords[0] == 0)
			{
				y = 0;
				for (int i = 0; i < Nx; i++)
				{
					x = start[0] + i*dx;
					FORMAT(vOut, x, y, v_down[i], 0.0, 0.0, 0.0);
				}
			}
			// Left: Located at the 'left-most' edge. Correspond to a coords[] column index = 0
			if (coords[1] == 0)
			{
				x = 0;
				for (int j = 0; j < Ny; j++)
				{
					y = start[1] + j*dy;
					FORMAT(vOut, x, y, v_left[j], 0.0, 0.0, 0.0);
				}
			}

			// Right: Located at the 'right-most' edge. Correspond to the last value in coords[] column index = number of partitions in x-direction
			if (coords[1] == Px - 1)
			{
				x = Lx;
				for (int j = 0; j < Ny; j++)
				{
					y = start[1] + j*dy;
					FORMAT(vOut, x, y, v_right[j], 0.0, 0.0, 0.0);
				}
			}
			vOut.close();
		}
		MPI_Barrier(mygrid);											//Blocks MPI caller until all processes within the mygrid communicator have reach this point of the code
	}
	if (rank == 0){
			cout<<"Data output terminated successfully" <<endl;
	}
}