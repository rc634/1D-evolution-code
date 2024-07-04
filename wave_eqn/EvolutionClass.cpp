#ifndef EVOLUTIONCLASS_CPP
#define EVOLUTIONCLASS_CPP

EvolutionClass::EvolutionClass(ParamsClass &a_params) :
           m_gridpoints(a_params.m_gridpoints), 
           m_x_min(a_params.m_x_min) , 
           m_x_max(a_params.m_x_max), 
           m_dx(a_params.m_dx),
           m_alpha(a_params.m_alpha),
           m_order(a_params.m_order),
           m_params(a_params),
           m_old(a_params),
           m_1(a_params),
           m_2(a_params),
           m_3(a_params),
           m_4(a_params),
           m_new(a_params),
           m_save_freq(a_params.m_save_freq)
{
    m_dt = m_alpha * m_dx;
    m_num_fields = m_old.m_fields.size();
    m_ghosts = m_order/2;

    // gaussian for psi, not pi (pi is initial momentum)
    m_old.m_phi.initialiseGaussian(0.,0.5,1.);
    m_old.m_psi.initialiseDGaussianDX(0.,0.5,1.);
    
}


void EvolutionClass::saveData()
{
    m_old.saveData();
}

void EvolutionClass::applyBC()
{
    for (int i=0; i<m_num_fields; i++)
    {
        m_old.m_fields[i]->applyBC_symmetric_4th();
        m_1.m_fields[i]->applyBC_symmetric_4th();
        m_2.m_fields[i]->applyBC_symmetric_4th();
        m_3.m_fields[i]->applyBC_symmetric_4th();
        m_4.m_fields[i]->applyBC_symmetric_4th();
        m_new.m_fields[i]->applyBC_symmetric_4th();

        // m_old.m_fields[i]->applyBC_periodic_4th();
        // m_1.m_fields[i]->applyBC_periodic_4th();
        // m_2.m_fields[i]->applyBC_periodic_4th();
        // m_3.m_fields[i]->applyBC_periodic_4th();
        // m_4.m_fields[i]->applyBC_periodic_4th();
        // m_new.m_fields[i]->applyBC_periodic_4th();

        // m_old.m_fields[i]->applyBC_jank();
        // m_1.m_fields[i]->applyBC_jank();
        // m_2.m_fields[i]->applyBC_jank();
        // m_3.m_fields[i]->applyBC_jank();
        // m_4.m_fields[i]->applyBC_jank();
        // m_new.m_fields[i]->applyBC_jank();

        // m_old.m_fields[i]->applyBC_zero();
        // m_1.m_fields[i]->applyBC_zero();
        // m_2.m_fields[i]->applyBC_zero();
        // m_3.m_fields[i]->applyBC_zero();
        // m_4.m_fields[i]->applyBC_zero();
        // m_new.m_fields[i]->applyBC_zero();
    }
}

// perform the needed amount of rk4 timesteps
void EvolutionClass::do_rk4_step(int timesteps)
{
    // dt to be used in rk4
    double h = m_dt;
    
    // number of timesteps
    for (int t=0; t<timesteps; t++)
    {
        //gooooo
        
        if (t%m_save_freq==0)
        {
            std::cout << "At time : " << t*h << std::endl;
            saveData();
        }        

        apply_RHS(m_old, m_4, m_1, 0.0); //m_4 is multiplied by 0 internally here
        applyBC();
        apply_RHS(m_old, m_1, m_2, 0.5 * h);
        applyBC();
        apply_RHS(m_old, m_2, m_3, 0.5 * h);
        applyBC();
        apply_RHS(m_old, m_3, m_4, h);
        applyBC();

        for (int n=0; n<m_num_fields; n++)
        {
            for (int i=m_ghosts; i<m_gridpoints-m_ghosts; i++)
            {
                m_new.m_fields[n]->m_field[i] = m_old.m_fields[n]->m_field[i];

                m_new.m_fields[n]->m_field[i] += h*( m_1.m_fields[n]->m_field[i] 
                                                      + m_4.m_fields[n]->m_field[i])/6.; 
                m_new.m_fields[n]->m_field[i] += h*( m_2.m_fields[n]->m_field[i] 
                                                      + m_3.m_fields[n]->m_field[i])/3.;
            
                m_old.m_fields[n]->m_field[i] = m_new.m_fields[n]->m_field[i];

                m_time += h;
            }
        }
        applyBC();
    }
}

// application of differential equation
void EvolutionClass::apply_RHS(PhysicsClass &old, 
                               PhysicsClass &input, 
                               PhysicsClass &output, 
                                      double delta)
{
    // simple wave equation
    for (int i=m_ghosts; i<m_gridpoints-m_ghosts; i++)
    {
        output.m_phi.m_field[i] = old.m_psi.m_field[i] + delta * input.m_psi.m_field[i];
        output.m_psi.m_field[i] = old.m_phi.d2_4th(i) + delta * input.m_phi.d2_4th(i);
    }
}


#endif // EVOLUTIONCLASS_CPP