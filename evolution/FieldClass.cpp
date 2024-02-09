#ifndef FIELDCLASS_CPP
#define FIELDCLASS_CPP

FieldClass::FieldClass(ParamsClass &a_params, std::string a_name) :
          m_gridpoints(a_params.m_gridpoints), 
          m_x_min(a_params.m_x_min), 
          m_x_max(a_params.m_x_max), 
          m_dx(a_params.m_dx),
          m_field_name(a_name),
          m_bc_type(a_params.m_bc_type),
          m_bc_pos(a_params.m_bc_pos),
          m_order(a_params.m_order)
{
    // Resize the vector to the specified length, and zero it
    m_field.resize(a_params.m_gridpoints,0.);
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
        double x = getX(i);
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
            double value = getValue(i);
            outputFile << getX(i) << " " << value << std::endl;
        }
        outputFile.close();
        //std::cout << "DAT file saved successfully." << std::endl;
    } 
    else 
    {
        std::cerr << "Error: Unable to save " + m_field_name + ".dat" << std::endl;
    }
}

// symmetric BC's for 4th order
void FieldClass::applyBC_symmetric_4th()
{
    int imax = m_gridpoints-1;
    if (m_bc_pos==0)
    {
        m_field[0] = m_field[3];
        m_field[1] = m_field[2];
        m_field[imax] = m_field[imax-3];
        m_field[imax-1] = m_field[imax-2];
    }
    else if (m_bc_pos==1)
    {
        m_field[0] = m_field[4];
        m_field[1] = m_field[3];
        m_field[imax] = m_field[imax-4];
        m_field[imax-1] = m_field[imax-3];
    }
}

double FieldClass::getX(int i) const
{
    if (m_bc_pos==0)
    {
        return (m_x_min - 0.5 * ((double) m_order - 1.) * m_dx) + i * m_dx;
    }
    else if (m_bc_pos==1)
    {
        return (m_x_min - 0.5 * m_order * m_dx) + i * m_dx;
    }
    else 
    {
        std::cout << "Improper m_bc_pos use in getX " << std::endl;
        return -7.7;
    }
}

#endif // FIELDCLASS_CPP