#ifndef FIELDCLASS_HPP
#define FIELDCLASS_HPP

#include <vector>
#include "ParamsClass.hpp"

class FieldClass {

public:

    std::vector<double> m_field;

private:
    int m_gridpoints;
    double m_x_min;
    double m_x_max;
    double m_dx, m_dx_sqr;
    int m_bc_type;
    int m_bc_pos;
    int m_order;
    std::string m_field_name;

    double m_twelfth = 1./12.;
    double m_twothirds = 2./3.;
    double m_fourthirds = 4./3.;

public:
    // Constructor that initializes the vector with a specified length
    FieldClass(ParamsClass &a_params, std::string a_name);

    // Function to set a value at a specific index in the vector
    void setValue(int index, double value);

    // Function to get the value at a specific index in the vector
    double getValue(int index) const;

    // initialises array with gaussian
    void initialiseGaussian(double mean, double sigma, double a);

    // initialises array with d/dx gaussian(x) 
    void initialiseDGaussianDX(double mean, double sigma, double a);

    // save data to dat file 
    void saveData() const;

    // apply symmetric boundary conditions
    void applyBC_symmetric_4th();

    void applyBC_zero();

    void applyBC_periodic_4th();

    void applyBC_jank();

    // fetch x coordinate at point i, taking into acount ghosts adn boundary positions
    double getX(int i) const; 

    // 4th order 1st derivative
    inline double d1_4th(int i);

    // 4th order 2nd derivative
    inline double d2_4th(int i);

    // 2nd order 2nd derivative
    inline double d2_2nd(int i);

};


#include "FieldClass.cpp"

#endif // FIELDCLASS_HPP