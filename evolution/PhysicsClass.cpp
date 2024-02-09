#ifndef PHYSICSCLASS_CPP
#define PHYSICSCLASS_CPP

PhysicsClass::PhysicsClass(int length, double a_x_min, double a_x_max,  
                                  std::vector<FieldClass> &a_all_fields) :
           m_gridpoints(length), 
           m_x_min(a_x_min) , 
           m_x_max(a_x_max), 
           m_all_fields(a_all_fields),
           m_psi(a_all_fields[0]),
           m_pi(a_all_fields[1])
{
    m_dx = (m_x_max-m_x_min)/((double) m_gridpoints-1);
}




void PhysicsClass::saveData()
{
    for (const auto field_obj : m_all_fields )
    {
        field_obj.saveData();
    }
}

#endif // PHYSICSCLASS_CPP