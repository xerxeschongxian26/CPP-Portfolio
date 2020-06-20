#pragma once //Preprocessor directive so current source file is included only once in a single compilation;improve compilation speed	
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options; //Namespace instance defined to reduce unnecessary repetition

/*Declarations for functions found in Prog_Options.cpp;
============================================================================================================================================================================ */

bool OptionStatus(int argc, char*argv[], po::variables_map &vm);

void ParamReader(po::variables_map &vm, double &Lx, double &Ly, int &Nx, int &Ny, int &Px, int &Py, double &dt, double &T, double &Re);