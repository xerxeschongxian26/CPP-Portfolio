#include <iostream>
#include <exception>
#include <iomanip>
#include <cmath>

#include "LidDrivenCavity.h"
#include "Prog_Options.h"
#include "MPI_Config.h"
#include <mpi.h>

using namespace std;

int main(int argc, char **argv)
{
/*MPI Initialise and Check;
================================================================================================================================ */

    int rank=0;
    int size=0;
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        cout << "MPI Initialisation Failed" << endl;
        return -1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 						// Obtain rank AKA proc ID
    MPI_Comm_size(MPI_COMM_WORLD, &size); 						// Obtain total number of processors.

/*Parsing of command line arguments;
================================================================================================================================== */

    //Variables for parsing of command line arguments using boost/program_options
    po::variables_map vm;      
    bool status_op,status_proc,status_cfl;
    
    //Variables for MPI Implementation
    int Neighbour[4] = {0};                                                            // 4 Neighbours elements initialised at 0
    int nx; 									       //Size of local x-grid for each rank
    int ny; 									       //Size of local y-grid for each rank										
    double start[2] = {0.0, 0.0};

    //Predefined variables for command line arguments
    double Lx; 
    double Ly;
    int    Nx;
    int    Ny;
    int    Px;
    int    Py;
    double dt;
    double T;
    double Re;

/*Checks for program options, correct num of processes, CFL condition for time step;
================================================================================================================================= */

    status_op = OptionStatus(argc, argv, vm); 						//Return 0 if fail, 1 if successfull
    if (status_op==0) {
        if (rank == 0) {
            cout << "An error has occurred in running the program options library" << endl;
        }
        MPI_Finalize();
        return 0;
    }
    ParamReader(vm, Lx, Ly, Nx, Ny, Px, Py, dt, T, Re);
    double dx = Lx/(Nx - 1);
    double dy = Ly/(Ny - 1);

    status_proc=ProcessCheck(size, Px, Py);					   
    if (status_proc==0) {
        if (rank == 0) {
            cout << "Number of processors not equal Px*Py." << endl;;
        }
        MPI_Finalize();
        return 0;
    }

    status_cfl=cflCheck(dt, Re, dx, dy);
    if (status_cfl==0) {
        if (rank == 0) {
           cout << "Timestep does not satisfy CFL condition. Please use different inputs" << endl;
        }
        MPI_Finalize();
        return 0;
    }
   
/*Initialise Cartesian Grid Topology;
======================================================================================================================================== */

    MPI_Comm mygrid;									//New communicator defined
    const int dims    = 2;  
    int sizes[dims]   = {Py, Px}; 							//2D grid with 2 coordinates: y(row) and x(col)
    int periods[dims] = {0, 0}; 							//Non periodic grid
    int reorder       = 0;								//No reordering 
    
    // New communicator from Cartesian Topology
    MPI_Cart_create(MPI_COMM_WORLD, dims, sizes, periods, reorder, &mygrid);
    int coords[dims];

    /// Compute coordinates of this process rank and store in coords
    MPI_Cart_coords(mygrid, rank, dims, coords);

    /// Obtain neighbour at each rank
    GetNeighbour(mygrid, Neighbour);

    //Staffer distributes work to ranks
    Staffer(rank, Nx, Ny, Px, Py, dx, dy, coords, start, nx, ny);    

/*Solver Implementation;
======================================================================================================================================== */

    LidDrivenCavity* solver = new LidDrivenCavity(mygrid, rank, coords, start, Neighbour,nx, ny, dt, T, Re, dx, dy);
    solver->FULLSOLVER();

/*Output to file and Termination;
======================================================================================================================================== */
    solver->outputFile(Lx, Ly, Px, Py);

    MPI_Finalize();
    return 0;
}




