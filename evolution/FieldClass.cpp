#ifndef FIELDCLASS_CPP
#define FIELDCLASS_CPP

FieldClass::FieldClass(ParamsClass &a_params, std::string a_name) :
          m_gridpoints(a_params.m_gridpoints), 
          m_x_min(a_params.m_x_min), 
          m_x_max(a_params.m_x_max), 
          m_dx(a_params.m_dx),
          m_dx_sqr(a_params.m_dx*a_params.m_dx),
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

void FieldClass::initialiseDGaussianDX(double mean, double sigma, double a)
{
    for (int i = 0; i < m_gridpoints; i++) 
    {
        double x = getX(i);
        m_field[i] = a*exp(-0.5 * pow((x - mean) / sigma, 2.0));

        // extra derivative for the velocity
        m_field[i] *= (x - mean)/(sigma*sigma);
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

// void FieldClass::saveData() const 
// {
//     std::ofstream outputFile(m_field_name + ".dat");
//     if (outputFile.is_open()) 
//     {
//         // Write the data to the CSV file
//         for (int i = 0; i < m_gridpoints; i++) 
//         {
//             double value = getValue(i);
//             outputFile << getX(i) << " " << value << std::endl;
//         }
//         outputFile.close();
//         //std::cout << "DAT file saved successfully." << std::endl;
//     } 
//     else 
//     {
//         std::cerr << "Error: Unable to save " + m_field_name + ".dat" << std::endl;
//     }
// }

void FieldClass::saveData() const 
{
    std::ofstream outputFile(m_field_name + ".dat", std::ios_base::app);
    if (outputFile.is_open()) 
    {
        // Write the data to the CSV file
        for (int i = 0; i < m_gridpoints-1; i++) 
        {
            outputFile << m_field[i] << ",";
        }
        outputFile << m_field[m_gridpoints-1] << std::endl;
        outputFile.close();
        //std::cout << "DAT file saved successfully." << std::endl;
    } 
    else 
    {
        std::cerr << "Error: Unable to save " + m_field_name + ".dat" << std::endl;
    }
}

inline double FieldClass::d1_4th(int i)
{
    double deriv=0.;
    deriv += m_twelfth*(m_field[i-2] - m_field[i+2]);
    deriv += m_twothirds*(-m_field[i-1] + m_field[i+1]);
    return deriv/m_dx;
}

inline double FieldClass::d2_4th(int i)
{
    double deriv=0.;
    deriv += -m_twelfth*(m_field[i-2] + m_field[i+2]);
    deriv += m_fourthirds*(m_field[i-1] + m_field[i+1]);
    deriv += -2.5*m_field[i];
    return deriv/m_dx_sqr;
}

inline double FieldClass::d2_2nd(int i)
{
    double deriv=0.;
    deriv += m_field[i-1] - m_field[i+1];
    deriv += -2.*m_field[i];
    return deriv/(m_dx*m_dx);
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

// symmetric BC's for 4th order
void FieldClass::applyBC_periodic_4th()
{
    int imax = m_gridpoints-1;
    if (m_bc_pos==0)
    {
        m_field[0] = m_field[imax-3];
        m_field[1] = m_field[imax-2];
        m_field[imax] = m_field[3];
        m_field[imax-1] = m_field[2];
    }
    else if (m_bc_pos==1)
    {
        m_field[0] = m_field[imax-4];
        m_field[1] = m_field[imax-3];
        m_field[imax] = m_field[4];
        m_field[imax-1] = m_field[3];
    }
}

// test cursed boundary conditions
void FieldClass::applyBC_jank()
{
    int imax = m_gridpoints-1;
    if (m_bc_pos==0)
    {
        m_field[0] = m_field[3];
        m_field[1] = m_field[2];
        m_field[imax] = -7.;
        m_field[imax-1] = -8.;
    }
    else if (m_bc_pos==1)
    {
        m_field[0] = m_field[4];
        m_field[1] = m_field[3];
        m_field[imax] = m_field[imax-4];
        m_field[imax-1] = m_field[imax-3];
    }
}


void FieldClass::applyBC_zero()
{
    int imax = m_gridpoints-1;
    
    m_field[0] = 0.;
    m_field[1] = 0.;
    m_field[imax] = 0.;
    m_field[imax-1] = 0.;
    
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