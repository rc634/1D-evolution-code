#ifndef FIELDCLASS_CPP
#define FIELDCLASS_CPP

FieldClass::FieldClass(int length, double a_x_min, double a_x_max, std::string a_name) 
          : m_gridpoints(length), m_x_min(a_x_min) , m_x_max(a_x_max), m_field_name(a_name)
{
    // Resize the vector to the specified length, and zero it
    m_field.resize(length,0.);
    m_dx = (m_x_max-m_x_min)/((double) m_gridpoints-1);
}

void FieldClass::setValue(int index, double value) 
{
    // Set the value at the specified index in the vector
    if (index >= 0 && index < m_field.size()) 
    {
        m_field[index] = value;
    } 
    else 
    {
        // Handle index out of range error
        throw std::out_of_range("Index out of range");
    }
}

double FieldClass::getValue(int index) const
{
    // Get the value at the specified index in the vector
    if (index >= 0 && index < m_field.size()) 
    {
        return m_field[index];
    } 
    else 
    {
        // Handle index out of range error
        throw std::out_of_range("Index out of range");
    }
}

void FieldClass::initialiseGaussian(double mean, double sigma, double a)
{
    for (int i = 0; i < m_gridpoints; i++) 
    {
        double x = m_x_min + i * m_dx;
        m_field[i] = a*exp(-0.5 * pow((x - mean) / sigma, 2.0));
    }

}

void FieldClass::saveData() const 
{
    std::ofstream outputFile(m_field_name + ".dat");
    if (outputFile.is_open()) 
    {
        // Write the data to the CSV file
        for (int i = 0; i < m_gridpoints; i++) 
        {
            double x_ = m_x_min + i * m_dx;
            double value = getValue(i);
            outputFile << x_ << " " << value << std::endl;
        }
        outputFile.close();
        //std::cout << "DAT file saved successfully." << std::endl;
    } 
    else 
    {
        std::cerr << "Error: Unable to save " + m_field_name + ".dat" << std::endl;
    }
}

#endif // FIELDCLASS_CPP