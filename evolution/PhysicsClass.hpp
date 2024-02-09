#ifndef PHYSICSCLASS_HPP
#define PHYSICSCLASS_HPP

#include <vector>

class PhysicsClass {

private:

    
    int m_gridpoints;
    double m_x_min;
    double m_x_max;
    double m_dx;

    std::vector<FieldClass> &m_all_fields;
    FieldClass &m_psi;
    FieldClass &m_pi;

public:
    // Constructor that initializes the vector with a specified length
    PhysicsClass(int length, double a_x_min, double a_x_max, std::vector<FieldClass> &a_all_fields);

    // save data to dat file 
    void saveData();
};


#include "PhysicsClass.cpp"

#endif // PHYSICSCLASS_HPP