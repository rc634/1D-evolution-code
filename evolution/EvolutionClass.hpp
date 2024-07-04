#ifndef EVOLUTIONCLASS_HPP
#define EVOLUTIONCLASS_HPP

#include <vector>
#include "ParamsClass.hpp"

class EvolutionClass {

private:

    int m_gridpoints;
    double m_x_min;
    double m_x_max;
    double m_dx;
    double m_alpha;
    int m_order, m_ghosts; // ghosts is half order
    int m_num_fields;
    int m_save_freq;
    double m_dt;
    double m_time = 0.;
    ParamsClass m_params;

    // old is input set of all fields, initial data for each timestep
    PhysicsClass m_old;

    // runge kutte 4 states
    PhysicsClass m_1;
    PhysicsClass m_2;
    PhysicsClass m_3;
    PhysicsClass m_4;

    // new data, output of runge kutte step
    PhysicsClass m_new;

public:
    // Constructor that initializes the vector with a specified length
    EvolutionClass(ParamsClass &a_params);

    // save data to dat file 
    void saveData();

    // does a runge kutte 4 time step
    void do_rk4_step(int timesteps);

private:

    // application of differential equation
    void apply_RHS(PhysicsClass &old, 
                   PhysicsClass &input, 
                   PhysicsClass &ouput, 
                          double delta);

    void applyBC();

};


#include "EvolutionClass.cpp"

#endif // EVOLUTIONCLASS_HPP