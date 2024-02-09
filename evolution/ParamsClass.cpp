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
}

#endif // PARAMSCLASS_CPP