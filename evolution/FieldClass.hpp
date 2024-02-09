#ifndef FIELDCLASS_HPP
#define FIELDCLASS_HPP

#include <vector>

class FieldClass {

private:

    std::vector<double> m_field;
    int m_gridpoints;
    double m_x_min;
    double m_x_max;
    double m_dx;
    std::string m_field_name;

public:
    // Constructor that initializes the vector with a specified length
    FieldClass(int length, double a_x_min, double a_x_max, std::string a_name);

    // Function to set a value at a specific index in the vector
    void setValue(int index, double value);

    // Function to get the value at a specific index in the vector
    double getValue(int index) const;

    // initialises array with gaussian
    void initialiseGaussian(double mean, double sigma, double a);

    // save data to dat file 
    void saveData() const;
};


#include "FieldClass.cpp"

#endif // FIELDCLASS_HPP