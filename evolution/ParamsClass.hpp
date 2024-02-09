#ifndef PARAMSCLASS_HPP
#define PARAMSCLASS_HPP

#include <vector>

class ParamsClass {

public:

    // default params - do not rely on
    int m_gridpoints;
    double m_x_min;
    double m_x_max;
    std::string m_filename;

    // Constructor that initializes the vector with a specified length
    ParamsClass(std::string a_filename);

    // Function to set a value at a specific index in the vector
    void loadParams();
};


#include "ParamsClass.cpp"

#endif // PARAMSCLASS_HPP