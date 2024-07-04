#ifndef PARAMSCLASS_CPP
#define PARAMSCLASS_CPP

ParamsClass::ParamsClass(std::string a_filename)
                      : m_filename(a_filename)
{
}

void ParamsClass::loadParams() 
{
    std::ifstream paramsFile(m_filename);
    if (!paramsFile.is_open()) {
        std::cerr << "Error: Unable to open parameters file." << std::endl;
    }

    std::string line;

    while (std::getline(paramsFile, line)) 
    {
        std::istringstream iss(line);
        std::string paramName;
        if (iss >> paramName) 
        {
            if (paramName == "gridpoints") 
            {
                if (!(iss >> m_gridpoints)) 
                {
                    std::cerr << "Error: Invalid value for vector_gridpoints." << std::endl;
                }
            } 

            else if (paramName == "x_min") 
            {
                if (!(iss >> m_x_min)) 
                {
                    std::cerr << "Error: Invalid value for x_min." << std::endl;
                }
            } 
            
            else if (paramName == "x_max") 
            {
                if (!(iss >> m_x_max)) 
                {
                    std::cerr << "Error: Invalid value for x_max." << std::endl;
                }
            } 

            else if (paramName == "alpha") 
            {
                if (!(iss >> m_alpha)) 
                {
                    std::cerr << "Error: Invalid value for alpha." << std::endl;
                }
            } 

            else if (paramName == "order") 
            {
                if (!(iss >> m_order)) 
                {
                    std::cerr << "Error: Invalid value for m_order." << std::endl;
                }
            } 

            else if (paramName == "bc_position") 
            {
                if (!(iss >> m_bc_pos)) 
                {
                    std::cerr << "Error: Invalid value for BC position." << std::endl;
                }
            } 

            else if (paramName == "bc_type") 
            {
                if (!(iss >> m_bc_type)) 
                {
                    std::cerr << "Error: Invalid value for BC type." << std::endl;
                }
            } 

            else if (paramName == "n_timesteps") 
            {
                if (!(iss >> m_n_timesteps)) 
                {
                    std::cerr << "Error: Invalid value n_timesteps." << std::endl;
                }
            } 

            else if (paramName == "save_freq") 
            {
                if (!(iss >> m_save_freq)) 
                {
                    std::cerr << "Error: Invalid value save_freq." << std::endl;
                }
            } 
            
            else 
            {
                std::cerr << "Error: Unexpected parameter name in the parameters file." << std::endl;
            }
        }
    }
    paramsFile.close();

    // some error warnings - maybe expand this 

    if (m_gridpoints <= 0) 
    {
        std::cerr << "Error: Invalid value for gridpoints." << std::endl;
    }

    if (m_x_min >= m_x_max) 
    {
        std::cerr << "Error: x_min must be smaller than x_max." << std::endl;
    }

    // calculate determined parameters
    m_gridpoints_phys = m_gridpoints; // non-ghost gridpoints
    m_gridpoints += m_order;

    if (m_bc_pos == 0) // midpoint
    {
        m_dx = (m_x_max-m_x_min)/((double) m_gridpoints_phys);
    }
    else if (m_bc_pos == 1) //centred
    {
        m_dx = (m_x_max-m_x_min)/((double) m_gridpoints_phys -1 );
    }

    

}

#endif // PARAMSCLASS_CPP