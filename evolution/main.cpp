#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

// include field class for physical data
#include "FieldClass.hpp" 
#include "ParamsClass.hpp"
#include "PhysicsClass.hpp"
#include "EvolutionClass.hpp"


//clean working directory of old data
void deleteLocalDatFiles() {
    std::string directory = "./"; // Specify the directory here

    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.path().extension() == ".dat") {
            std::filesystem::remove(entry.path());
            std::cout << "Deleted: " << entry.path() << std::endl;
        }
    }
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

    //////////////
    // delete dat files in dir

    deleteLocalDatFiles();


    ///////////////
    // params 

    ParamsClass params(argv[1]);
    params.loadParams();

    /////////////
    // evolve PDE
    
    EvolutionClass waveEQN(params);

    waveEQN.do_rk4_step(params.m_n_timesteps);

    ////////////
    // data saved automatically


    ///////////
    // end
    
    return 0;
}