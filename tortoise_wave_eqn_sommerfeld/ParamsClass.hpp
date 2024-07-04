#ifndef PARAMSCLASS_HPP
#define PARAMSCLASS_HPP

#include <vector>

class ParamsClass {

public:

    // default params - do not rely on
    int m_gridpoints;
    int m_gridpoints_phys;
    double m_x_min;
    double m_x_max;
    double m_alpha;
    std::string m_filename;
    double m_dx;
    int m_order;
    int m_bc_pos;
    int m_bc_type;
    int m_n_timesteps;
    int m_save_freq;
    double m_r_ext;

    // Constructor that initializes the vector with a specified length
    ParamsClass(std::string a_filename);

    // Function to set a value at a specific index in the vector
    void loadParams();
};


#include "ParamsClass.cpp"

#endif // PARAMSCLASS_HPP
