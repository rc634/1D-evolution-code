#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

// include field class for physical data
#include "FieldClass.hpp" 


// Temporary gaussian function here 
// Function to calculate Gaussian wavepacket
double gaussianWavepacket(double x, double mean, double sigma) {
    return exp(-0.5 * pow((x - mean) / sigma, 2.0));
}





int main() 
{
    ///////////////
    // params 

    std::ifstream paramsFile("params.txt");
    if (!paramsFile.is_open()) {
        std::cerr << "Error: Unable to open parameters file." << std::endl;
        return 1;
    }

    // data to be read from params.txt
    int gridpoints = 0;
    double x_min = 0.0, x_max = 0.0;
    std::string line;

    while (std::getline(paramsFile, line)) 
    {
        std::istringstream iss(line);
        std::string paramName;
        if (iss >> paramName) 
        {
            if (paramName == "gridpoints") 
            {
                if (!(iss >> gridpoints)) 
                {
                    std::cerr << "Error: Invalid value for vector_gridpoints." << std::endl;
                    return 1;
                }
            } else if (paramName == "x_min") 
            {
                if (!(iss >> x_min)) 
                {
                    std::cerr << "Error: Invalid value for x_min." << std::endl;
                    return 1;
                }
            } 
            else if (paramName == "x_max") 
            {
                if (!(iss >> x_max)) 
                {
                    std::cerr << "Error: Invalid value for x_max." << std::endl;
                    return 1;
                }
            } 
            else 
            {
                std::cerr << "Error: Unexpected parameter name in the parameters file." << std::endl;
                return 1;
            }
        }
    }
    paramsFile.close();

    if (gridpoints <= 0) {
        std::cerr << "Error: Invalid value for gridpoints." << std::endl;
        return 1;
    }
    

    ///////////
    //init

    const double dx = (x_max - x_min) / (gridpoints - 1); // Step size

    FieldClass psi(gridpoints);
    
    ////////////////
    // Initialize psi with Gaussian wavepacket

    double mean = 0.0; // Mean of the Gaussian
    double sigma = 1.0; // Standard deviation of the Gaussian
    for (int i = 0; i < gridpoints; ++i) 
    {
        double x = x_min + i * dx;
        double value = gaussianWavepacket(x, mean, sigma);
        psi.setValue(i, value);
    }

    /////////////
    // Save the array to a CSV DAT file

    std::ofstream outputFile("wavepacket.dat");
    if (outputFile.is_open()) 
    {
        // Write the data to the CSV file
        for (int i = 0; i < gridpoints; ++i) 
        {
            double x_ = x_min + i * dx;
            double value = psi.getValue(i);
            outputFile << x_ << " " << value << std::endl;
        }
        outputFile.close();
        std::cout << "DAT file saved successfully." << std::endl;
    } 
    else 
    {
        std::cerr << "Error: Unable to open the DAT file." << std::endl;
        return 1;
    }

    return 0;
}