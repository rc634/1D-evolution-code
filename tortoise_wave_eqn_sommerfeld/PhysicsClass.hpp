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
    int m_num_fields;
    int m_bc_pos;
    double m_dt;
    double m_time = 0.;
    ParamsClass m_params;

public:

    std::vector<FieldClass*> m_fields;

    FieldClass m_phi;
    FieldClass m_psi;





    // Constructor that initializes the vector with a specified length
    PhysicsClass(ParamsClass &a_params);

    // save data to dat file
    void saveData();

};


#include "PhysicsClass.cpp"

#endif // PHYSICSCLASS_HPP
