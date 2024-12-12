#ifndef PHYSICSCLASS_CPP
#define PHYSICSCLASS_CPP

PhysicsClass::PhysicsClass(ParamsClass &a_params) :
           m_gridpoints(a_params.m_gridpoints),
           m_x_min(a_params.m_x_min) ,
           m_x_max(a_params.m_x_max),
           m_dx(a_params.m_dx),
           m_alpha(a_params.m_alpha),
           m_order(a_params.m_order),
           m_params(a_params),
           m_U(a_params,"U"),
           m_dUdt(a_params,"dUdt"),
           m_H(a_params,"H"),
           m_dHdt(a_params,"dHdt"),
           m_bc_pos(a_params.m_bc_pos)
{
    m_dt = m_alpha * m_dx;

    m_fields.push_back(&m_U);
    m_fields.push_back(&m_dUdt);
    m_fields.push_back(&m_H);
    m_fields.push_back(&m_dHdt);

    // DO NOT ADD X (coordinate) HERE AS IT WILL BE ITERATED THROUGH IN LATER

    m_num_fields = m_fields.size();
}


void PhysicsClass::saveData()
{
    for (const auto field_obj : m_fields )
    {
        field_obj->saveData();
    }
}


#endif // PHYSICSCLASS_CPP
