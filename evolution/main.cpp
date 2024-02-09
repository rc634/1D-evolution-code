#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

// include field class for physical data
#include "FieldClass.hpp" 
#include "ParamsClass.hpp"
#include "PhysicsClass.hpp"


// Temporary gaussian function here 
// Function to calculate Gaussian wavepacket
double gaussianWavepacket(double x, double mean, double sigma) 
{
    return exp(-0.5 * pow((x - mean) / sigma, 2.0));
}





int main(int argc, char* argv[]) 
{
    // passes name of parameter file ot argv[1]
    if (argc != 2) // allows only one argument to be passed
    {
        std::cerr << "Usage: " << argv[0] << " filename" << std::endl;
        return 1;
    }
    std::cout << "Parameter file name : " << argv[1] << std::endl;


    ///////////////
    // params 

    ParamsClass params(argv[1]);
    params.loadParams();
    

    ///////////
    // initialise each field


    // initialise fields
    FieldClass psi(params, "psi");
    FieldClass pi(params, "pi");

    // gaussian for psi, not pi (pi is initial momentum)
    psi.initialiseGaussian(0.,1.,3.);

    // set of all fields
    std::vector<FieldClass> all_fields;

    // add fields to set of all fields
    all_fields.push_back(psi);
    all_fields.push_back(pi);
    
    // add all fields to physics class
    PhysicsClass waveEQN(params, all_fields);

    //save data
    waveEQN.saveData();
    



    return 0;
}