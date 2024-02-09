#ifndef FIELDCLASS_CPP
#define FIELDCLASS_CPP

FieldClass::FieldClass(int length) {
    // Resize the vector to the specified length
    field.resize(length);
}

void FieldClass::setValue(int index, double value) {
    // Set the value at the specified index in the vector
    if (index >= 0 && index < field.size()) {
        field[index] = value;
    } else {
        // Handle index out of range error
        throw std::out_of_range("Index out of range");
    }
}

double FieldClass::getValue(int index) {
    // Get the value at the specified index in the vector
    if (index >= 0 && index < field.size()) {
        return field[index];
    } else {
        // Handle index out of range error
        throw std::out_of_range("Index out of range");
    }
}

#endif // FIELDCLASS_CPP