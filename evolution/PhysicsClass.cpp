#ifndef PHYSICSCLASS_CPP
#define PHYSICSCLASS_CPP

PhysicsClass::PhysicsClass(ParamsClass &a_params,  
                                  std::vector<FieldClass> &a_all_fields) :
           m_gridpoints(a_params.m_gridpoints), 
           m_x_min(a_params.m_x_min) , 
           m_x_max(a_params.m_x_max), 
           m_dx(a_params.m_dx),
           m_alpha(a_params.m_alpha),
           m_all_fields(a_all_fields),
           m_psi(a_all_fields[0]),
           m_pi(a_all_fields[1]),
           m_order(a_params.m_order)
{
}


void PhysicsClass::saveData()
{
    for (const auto field_obj : m_all_fields )
    {
        field_obj.saveData();
    }
}

#endif // PHYSICSCLASS_CPP