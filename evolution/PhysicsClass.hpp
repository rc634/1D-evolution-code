#ifndef PHYSICSCLASS_HPP
#define PHYSICSCLASS_HPP

#include <vector>
#include "ParamsClass.hpp"

class PhysicsClass {

private:

    int m_gridpoints;
    double m_x_min;
    double m_x_max;
    double m_dx;
    double m_alpha;
    int m_order;

    std::vector<FieldClass> &m_all_fields;
    FieldClass &m_psi;
    FieldClass &m_pi;

public:
    // Constructor that initializes the vector with a specified length
    PhysicsClass(ParamsClass &a_params, std::vector<FieldClass> &a_all_fields);

    // save data to dat file 
    void saveData();

};


#include "PhysicsClass.cpp"

#endif // PHYSICSCLASS_HPP