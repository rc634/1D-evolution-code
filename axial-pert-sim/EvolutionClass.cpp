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
           m_save_freq(a_params.m_save_freq),
           m_x(a_params,"x"),
           m_V_UU(a_params,"V_UU"),
           m_V_UH(a_params,"V_UH"),
           m_V_HH(a_params,"V_HH"),
           m_bc_pos(a_params.m_bc_pos),
           m_bc_type(a_params.m_bc_type),
           m_r_ext(a_params.m_r_ext)
{
    m_dt = m_alpha * m_dx;
    m_num_fields = m_old.m_fields.size();
    m_ghosts = m_order/2;

    // gaussian for psi, not pi (pi is initial momentum)
    m_old.m_H.initialiseGaussian(600.,4.,1.,-1.);
    m_old.m_U.initialiseGaussian(600.,4.,0.1,-1.);
    m_old.m_dHdt.initialiseDGaussianDX(600.,4.,1.,-1.);
    m_old.m_dUdt.initialiseDGaussianDX(600.,4.,0.1,-1.);

    //algebraic initial conditions (needed for X)
    setXandV();

    //read numerical initial conditions (overrides V)
    readXandV();
}

void EvolutionClass::setXandV()
{
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
        m_V_UU.m_field[i] = 0.;
        m_V_UH.m_field[i] = 0.;
        m_V_HH.m_field[i] = 0.;
        //m_V.m_field[i] = x*x;
    }
}


void EvolutionClass::readXandV()
{
    // Open the TSV file
    //std::ifstream x_file("../initial_data/EMDBH_1/tort.dat");

    //std::string base = "a05_Q10073NEW";
    std::string path = "../initial_data/UH2/";
    std::string base = "a06_Q10187NEW";

    std::ifstream x_file(path+"tort_"+base+".dat");
    if (!x_file.is_open())
    {
        std::cerr << "Error opening x_file." << std::endl;
    }

    std::ifstream V_UU_file(path+"VUU_"+base+".dat");
    if (!V_UU_file.is_open())
    {
        std::cerr << "Error opening V_UU_file." << std::endl;
    }
    std::ifstream V_UH_file(path+"VUH_"+base+".dat");
    if (!V_UH_file.is_open())
    {
        std::cerr << "Error opening V_UH_file." << std::endl;
    }
    std::ifstream V_HH_file(path+"VHH_"+base+".dat");
    if (!V_HH_file.is_open())
    {
        std::cerr << "Error opening V_HH_file." << std::endl;
    }


    // Define a vector to store each dat file
    std::vector<double> V_UU_input;
    std::vector<double> V_UH_input;
    std::vector<double> V_HH_input;
    std::vector<double> x_input;

    // Read the x_file line by line
    std::string x_line;
    while (getline(x_file, x_line))
    {
        // Split the line into fields based on tab ('\t') delimiter
        std::stringstream ss(x_line);
        std::string field;
        while (getline(ss, field, '\t'))
        {
            // Store each field in the vector
            x_input.push_back(std::stod(field));
        }
    }

    // Close the x_file
    x_file.close();

    // Read the V_file line by line
    std::string V_line;
    while (getline(V_UU_file, V_line))
    {
        // Split the line into fields based on tab ('\t') delimiter
        std::stringstream ss(V_line);
        std::string field;
        while (getline(ss, field, '\t'))
        {
            // Store each field in the vector
            V_UU_input.push_back(std::stod(field));
        }
    }
    while (getline(V_UH_file, V_line))
    {
        // Split the line into fields based on tab ('\t') delimiter
        std::stringstream ss(V_line);
        std::string field;
        while (getline(ss, field, '\t'))
        {
            // Store each field in the vector
            V_UH_input.push_back(std::stod(field));
        }
    }
    while (getline(V_HH_file, V_line))
    {
        // Split the line into fields based on tab ('\t') delimiter
        std::stringstream ss(V_line);
        std::string field;
        while (getline(ss, field, '\t'))
        {
            // Store each field in the vector
            V_HH_input.push_back(std::stod(field));
        }
    }

    // Close the V_XX_files
    V_UU_file.close();
    V_UH_file.close();
    V_HH_file.close();

    std::cout << "x_input[0] = " << x_input[0] << std::endl;
    std::cout << "x_input[final] = " << x_input[x_input.size()-1] << std::endl;

    if (m_x_min < x_input[0])
    {
        std::cout << "WARNING! x_min param too small!" << std::endl;
    }
    if (m_x_max > x_input[x_input.size()-1])
    {
        std::cout << "WARNING! x_max param too big!" << std::endl;
    }
    if (x_input.size() != V_UU_input.size())
    {
        std::cout << "WARNING! x and vUU input arrays of different length" << std::endl;
    }
    if (x_input.size() != V_UH_input.size())
    {
        std::cout << "WARNING! x and vUH input arrays of different length" << std::endl;
    }
    if (x_input.size() != V_HH_input.size())
    {
        std::cout << "WARNING! x and vHH input arrays of different length" << std::endl;
    }

    fourthOrderInterpXandV(x_input,V_UU_input, m_V_UU);
    fourthOrderInterpXandV(x_input,V_UH_input, m_V_UH);
    fourthOrderInterpXandV(x_input,V_HH_input, m_V_HH);
}


void EvolutionClass::fourthOrderInterpXandV(std::vector<double> &x_input, std::vector<double> &V_input, FieldClass &V_output)
{
    // length/size of input arrays
    int n_s = x_input.size();
    // dx of input x array
    double dx_s = (x_input[n_s-1]-x_input[0])/((double) n_s - 1.);

    double x_min_s = x_input[0];

    // input[i-1], input[i] input[i+1] input[i+2] for f and v
    double fx1, fx2, fx3, fx4;
    double fv1, fv2, fv3, fv4;

    // coefficient of polynomial, A + Bx + Cx^2/2 + Dx^3/6
    double Ax, Bx, Cx, Dx;
    double Av, Bv, Cv, Dv;

    // fraction above (iter + 1/2), i.e. midpoint
    double a;
    int iter; // integer gridpoint in x_input

    for (int i = 0; i < m_x.m_field.size(); i++)
    {
        iter = (int) floor((m_x.m_field[i]-x_min_s)/dx_s);
        a = (m_x.m_field[i]-x_min_s)/dx_s
                   - floor((m_x.m_field[i]-x_min_s)/dx_s)-0.5;

        //std::cout << ">*<>*< : " << iter << " : " << a << " : " << x_input[iter] << std::endl;

        fx1 = x_input[iter-1];
        fx2 = x_input[iter];
        fx3 = x_input[iter+1];
        fx4 = x_input[iter+2];
        fv1 = V_input[iter-1];
        fv2 = V_input[iter];
        fv3 = V_input[iter+1];
        fv4 = V_input[iter+2];

        if (iter<1)
        {
            std::cout << "m_x_min is slightly too small" << std::endl;
        }
        if (iter>x_input.size()-3)
        {
            std::cout << "m_x_max param is slightly too large" << std::endl;
        }

        Ax = (9.*(fx3+fx2)-(fx1+fx4))/16.;
        Bx = (27.*(fx3-fx2)-(fx4-fx1))/24.;
        Cx = (fx4+fx1-fx2-fx3)/2.;
        Dx = (fx4-fx1) - 3.*(fx3-fx2);
        Av = (9.*(fv3+fv2)-(fv1+fv4))/16.;
        Bv = (27.*(fv3-fv2)-(fv4-fv1))/24.;
        Cv = (fv4+fv1-fv2-fv3)/2.;
        Dv = (fv4-fv1) - 3.*(fv3-fv2);

        m_x.m_field[i] = Ax + Bx*a + Cx*a*a/2. + Dx*a*a*a/6.;
        V_output.m_field[i] = Av + a*(Bv + a*(Cv/2. + a*Dv/6.));
    }
}

// void EvolutionClass::secondOrderInterpXandV(std::vector<double> &x_input, std::vector<double> &V_input)
// {
//     // length/size of input arrays
//     int n_s = x_input.size();
//     // dx of input x array
//     double dx_s = (x_input[n_s-1]-x_input[0])/((double) n_s - 1.);
//
//     double x_min_s = x_input[0];
//
//     // input[i-1], input[i] input[i+1] input[i+2]
//     double fx1, fx2;
//     double fv1, fv2;
//
//     // coefficient of polynomial, A + Bx + Cx^2/2 + Dx^3/6
//     double Ax, Bx;
//     double Av, Bv;
//
//     // fraction above (iter + 1/2), i.e. midpoint
//     double a;
//     int iter; // integer gridpoint in x_input
//
//     for (int i = 0; i < m_x.m_field.size(); i++)
//     {
//         iter = (int) floor((m_x.m_field[i]-x_min_s)/dx_s);
//         a = (m_x.m_field[i]-x_min_s)/dx_s
//                    - floor((m_x.m_field[i]-x_min_s)/dx_s)-0.5;
//
//         std::cout << ">*<>*< : " << iter << " : " << a << " : " << x_input[iter] << std::endl;
//
//         fx1 = x_input[iter];
//         fx2 = x_input[iter+1];
//         fv1 = V_input[iter];
//         fv2 = V_input[iter+1];
//
//         if (iter<1)
//         {
//             std::cout << "m_x_min is slightly too small" << std::endl;
//         }
//         if (iter>x_input.size()-3)
//         {
//             std::cout << "m_x_max param is slightly too large" << std::endl;
//         }
//
//         Ax = (fx1+fx2)/2.;
//         Bx = (fx2-fx1);
//         Av = (fv1+fv2)/2.;
//         Bv = (fv2-fv1);
//
//         m_x.m_field[i] = Ax + Bx*a;
//         m_V.m_field[i] = Av + a*Bv;
//     }
// }


void EvolutionClass::saveData()
{
    m_old.saveData();
    m_x.saveData();
    m_V_UU.saveData();
    m_V_UH.saveData();
    m_V_HH.saveData();
}

void EvolutionClass::writeTimeDataPoint(std::string a_filename, double a_value)
{
    std::ofstream outputFile(a_filename, std::ios_base::app);
    if (outputFile.is_open())
    {
        outputFile << a_value << std::endl;
        outputFile.close();
        //std::cout << "DAT file saved successfully." << std::endl;
    }
}


double EvolutionClass::getInterpValue(std::vector<double> &a_field, double r)
{
    if (r + m_dx > m_x_max)
    {
        std::cout << "Warning, Interpolating above m_x_max!" << std::endl;
    }
    else if (r - m_dx < m_x_min)
    {
        std::cout << "Warning, Interpolating below m_x_min!" << std::endl;
    }

    // nearest integer place to r in m_x
    int iter = (int) floor((r-m_x.m_field[0])/m_dx);

    // the fractional integer corresponding to r,
    // the 0.5 means that in the middle of two values, delta_i=0.
    double delta_i = (r-m_x.m_field[0])/m_dx - ((double) iter) - 0.5;

    double f_m1 = a_field[iter-1];
    double f_0 = a_field[iter];
    double f_1 = a_field[iter+1];
    double f_2 = a_field[iter+2];

    double f_interp = 0.0625 * ( 9.*(f_0+f_1) - (f_2 + f_m1));
    f_interp += delta_i*(27.*(f_1-f_0) - (f_2-f_m1))/24.;
    f_interp += delta_i*delta_i*(f_2-f_1-f_0+f_m1)*0.25;
    f_interp += pow(delta_i,3)*(f_2-f_m1-3.*(f_1-f_0))/6.;

    return f_interp;
}


void EvolutionClass::saveTimeData()
{
    double _U_val = getInterpValue(m_old.m_U.m_field, m_r_ext);
    double _H_val = getInterpValue(m_old.m_H.m_field, m_r_ext);
    double _log_U = 0.5 * log(_U_val*_U_val);
    double _log_H = 0.5 * log(_H_val*_H_val);


    writeTimeDataPoint("U_of_time.dat", _U_val);
    writeTimeDataPoint("H_of_time.dat", _H_val);
    writeTimeDataPoint("log_U_of_time.dat", _log_U);
    writeTimeDataPoint("log_H_of_time.dat", _log_H);
    writeTimeDataPoint("time.dat", m_time);
}

void EvolutionClass::applyBC()
{

    // reflective/symmetric BCs
    if (m_bc_type==0)
    {
        for (int i=0; i<m_num_fields; i++)
        {
            m_old.m_fields[i]->applyBC_symmetric_4th();
            m_1.m_fields[i]->applyBC_symmetric_4th();
            m_2.m_fields[i]->applyBC_symmetric_4th();
            m_3.m_fields[i]->applyBC_symmetric_4th();
            m_4.m_fields[i]->applyBC_symmetric_4th();
            m_new.m_fields[i]->applyBC_symmetric_4th();
        }
        m_V_UU.applyBC_symmetric_4th();
        m_V_UH.applyBC_symmetric_4th();
        m_V_UU.applyBC_symmetric_4th();
    }

    // periodic BCs
    if (m_bc_type==1)
    {
        for (int i=0; i<m_num_fields; i++)
        {
            m_old.m_fields[i]->applyBC_periodic_4th();
            m_1.m_fields[i]->applyBC_periodic_4th();
            m_2.m_fields[i]->applyBC_periodic_4th();
            m_3.m_fields[i]->applyBC_periodic_4th();
            m_4.m_fields[i]->applyBC_periodic_4th();
            m_new.m_fields[i]->applyBC_periodic_4th();
        }
        m_V_UU.applyBC_periodic_4th();
        m_V_UH.applyBC_periodic_4th();
        m_V_HH.applyBC_periodic_4th();
    }

    // sommerfeld boundary conditions
    if (m_bc_type==2)
    {
            //pass, nothing to do here
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

        // save extraction radius data (lightweight)
        saveTimeData();

        // applyBC does nothing for Sommerfled, but is needed for other boundary types

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
            for (int i=0; i<m_gridpoints; i++)
            {
                m_new.m_fields[n]->m_field[i] = m_old.m_fields[n]->m_field[i];

                m_new.m_fields[n]->m_field[i] += h*( m_1.m_fields[n]->m_field[i]
                                                      + m_4.m_fields[n]->m_field[i])/6.;
                m_new.m_fields[n]->m_field[i] += h*( m_2.m_fields[n]->m_field[i]
                                                      + m_3.m_fields[n]->m_field[i])/3.;

                m_old.m_fields[n]->m_field[i] = m_new.m_fields[n]->m_field[i];
            }
        }
        applyBC();
        m_time += h;
    }
}

// application of differential equation
void EvolutionClass::apply_RHS(PhysicsClass &old,
                               PhysicsClass &input,
                               PhysicsClass &output,
                                      double delta)
{
    double LAMBDA = 1.; // damping params for ghosts



    // bulk evolution scheme
    for (int i=m_ghosts; i<m_gridpoints-m_ghosts; i++)
    {
        output.m_U.m_field[i] = old.m_dUdt.m_field[i] + delta * input.m_dUdt.m_field[i];
        output.m_H.m_field[i] = old.m_dHdt.m_field[i] + delta * input.m_dHdt.m_field[i];

        output.m_dUdt.m_field[i] = old.m_U.d2_4th(i) + delta * input.m_U.d2_4th(i);
        output.m_dUdt.m_field[i] += -m_V_UU.m_field[i] *
                                      (old.m_U.m_field[i] + delta * input.m_U.m_field[i]);
        output.m_dUdt.m_field[i] += -m_V_UH.m_field[i] *
                                      (old.m_H.m_field[i] + delta * input.m_H.m_field[i]);

        output.m_dHdt.m_field[i] = old.m_H.d2_4th(i) + delta * input.m_H.d2_4th(i);
        output.m_dHdt.m_field[i] += -m_V_UH.m_field[i] *
                                      (old.m_U.m_field[i] + delta * input.m_U.m_field[i]);
        output.m_dHdt.m_field[i] += -m_V_HH.m_field[i] *
                                      (old.m_H.m_field[i] + delta * input.m_H.m_field[i]);
    }

    if (m_bc_type==2) //sommerfeld = 2
    {
        // WAVE SPEED = 1 FOR SOMMERFELD BOUNDARY CONDITIONS

        // LHS boundary sommerfeld, ingoing wave d/dt - (1/v) d/dx = 0
        for (int i=0; i<m_ghosts; i++)
        {
            output.m_U.m_field[i] = old.m_U.d1_4th_right(i) + delta * input.m_U.d1_4th_right(i);
            output.m_H.m_field[i] = old.m_H.d1_4th_right(i) + delta * input.m_H.d1_4th_right(i);
            output.m_dUdt.m_field[i] = old.m_dUdt.d1_4th_right(i) + delta * input.m_dUdt.d1_4th_right(i);
            output.m_dHdt.m_field[i] = old.m_dHdt.d1_4th_right(i) + delta * input.m_dHdt.d1_4th_right(i);
        }

        // RHS boundary sommerfeld, outgoing wave d/dt + (1/v) d/dx = 0
        for (int i=m_gridpoints-m_ghosts; i<m_gridpoints; i++)
        {
          output.m_U.m_field[i] = -old.m_U.d1_4th_left(i) - delta * input.m_U.d1_4th_left(i);
          output.m_H.m_field[i] = -old.m_H.d1_4th_left(i) - delta * input.m_H.d1_4th_left(i);
          output.m_dUdt.m_field[i] = -old.m_dUdt.d1_4th_left(i) - delta * input.m_dUdt.d1_4th_left(i);
          output.m_dHdt.m_field[i] = -old.m_dHdt.d1_4th_left(i) - delta * input.m_dHdt.d1_4th_left(i);
        }
    }
}


#endif // EVOLUTIONCLASS_CPP
