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
           m_phi(a_params,"phi"),
           m_psi(a_params,"psi"),
           m_x(a_params,"x"),
           m_bc_pos(a_params.m_bc_pos)
{
    m_dt = m_alpha * m_dx;
    
    m_fields.push_back(&m_phi);
    m_fields.push_back(&m_psi);
    // DO NOT ADD X HERE AS IT WILL BE ITERATED THROUGH IN LATER

    m_num_fields = m_fields.size();

    double x, i_ng; // i_ng = i - Num_ghost
    for (int i=0; i<m_gridpoints; i++)
    {
        i_ng = ((double) i) - ((double) (m_order/2));
        if (m_bc_pos==0) // bc at midpoint 
        {
            x = m_x_min + (i_ng + 0.5) * m_dx;
        }
        else if (m_bc_pos==1) // bc at gridpoint
        {
            x = m_x_min + i_ng * m_dx;
        } 
        m_x.m_field[i] = x;
    }
}


void PhysicsClass::saveData()
{
    for (const auto field_obj : m_fields )
    {
        field_obj->saveData();
    }
    // must save x separately as it is not included in m_fields
    m_x.saveData();
}


void PhysicsClass::InitialiseArrays()
{

}




#endif // PHYSICSCLASS_CPP