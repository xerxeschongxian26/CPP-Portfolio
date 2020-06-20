#include<cstdlib>
#include<iostream>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;  //Class instance defined to reduce unnecessary repetition

/*Functions that control the parsing of command line arguments;
==================================================================================================================================== */

bool OptionStatus(int argc, char *argv[], po::variables_map &vm)
{
   try
   {
    po::options_description desc("Allowed input options");
    // Adding input options with appropriate default values
    desc.add_options()
        ("help", "Produce help message")
        ("Lx", po::value<double>()->default_value(1.0), "Length of domain in x-direction")
        ("Ly", po::value<double>()->default_value(1.0), "Length of domain in y-direction")
        ("Nx", po::value<int>()   ->default_value(161), "Number of grid points in x-direction")
        ("Ny", po::value<int>()   ->default_value(161), "Number of grid points in y-direction")
        ("Px", po::value<int>()   ->default_value(3), "Number of partitions in x-direction")
        ("Py", po::value<int>()   ->default_value(2), "Number of partitions in y-direction")
        ("dt", po::value<double>()->default_value(0.0005),"Time step size")
        ("T", po::value<double>() ->default_value(10000),"Final time")
        ("Re", po::value<double>()->default_value(100.0),"Reynolds")
    ;
    // Parse command-line arguments and store in buffer vm
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
   }
   catch(exception const &e)
   {
    cout << e.what() << endl;
    return false;
   }
   return true;
}

/*Functions that assigns the value from the command line arguments using the variable map vm created above;
====================================================================================================================================== */

void ParamReader(po::variables_map &vm, double &Lx, double &Ly, int &Nx, int &Ny, int &Px, int &Py, double &dt, double &T, double &Re)
{
    Lx = vm["Lx"].as<double>();
    Ly = vm["Ly"].as<double>();
    Nx = vm["Nx"].as<int>();
    Ny = vm["Ny"].as<int>();
    Px = vm["Px"].as<int>();
    Py = vm["Py"].as<int>();
    dt = vm["dt"].as<double>();
    T =  vm["T"] .as<double>();
    Re = vm["Re"].as<double>();
}