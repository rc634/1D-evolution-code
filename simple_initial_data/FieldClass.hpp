#ifndef FIELDCLASS_HPP
#define FIELDCLASS_HPP

#include <vector>

class FieldClass {

private:

    std::vector<double> field;

public:
    // Constructor that initializes the vector with a specified length
    FieldClass(int length);

    // Function to set a value at a specific index in the vector
    void setValue(int index, double value);

    // Function to get the value at a specific index in the vector
    double getValue(int index);
};


#include "FieldClass.cpp"

#endif // FIELDCLASS_HPP