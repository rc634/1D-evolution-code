#ifndef EVOLUTIONCLASS_HPP
#define EVOLUTIONCLASS_HPP

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
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
    int m_bc_pos;
    int m_bc_type;
    int m_save_freq;
    double m_dt;
    double m_time = 0.;
    ParamsClass m_params;
    double m_r_ext;

    // old is input set of all fields, initial data for each timestep
    PhysicsClass m_old;

    // runge kutte 4 states
    PhysicsClass m_1;
    PhysicsClass m_2;
    PhysicsClass m_3;
    PhysicsClass m_4;

    // new data, output of runge kutte step
    PhysicsClass m_new;


    // This is a field class to make use of the save functions
    // but it is not included into m_fields to avoid it's
    // automatic use later in the rk4 steps
    FieldClass m_x;
    FieldClass m_V;

public:
    // Constructor that initializes the vector with a specified length
    EvolutionClass(ParamsClass &a_params);

    // save data to dat file
    void saveData();

    // save data for phi(t,r=r_extraction) and pi
    void saveTimeData();

    // used in saveTimeData(), saves a single piece of data to a datafile
    void writeTimeDataPoint(std::string a_filename, double a_value);

    // 4th order interpolation to get field(x) at non-integer gridpoint
    double getInterpValue(std::vector<double> &a_field, double r);

    // does a runge kutte 4 time step
    void do_rk4_step(int timesteps);

    // setup array of x and V(x)
    void setXandV();
    void readXandV();
    void secondOrderInterpXandV(std::vector<double> &x_input, std::vector<double> &V_input);
    void fourthOrderInterpXandV(std::vector<double> &x_input, std::vector<double> &V_input);


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
