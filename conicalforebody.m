function co = conicalforebody(theta_c1, m_0, t_0, p_0, r_0,g_0)

        %function to calculate shock angle and the equivalent propertiesafter a conical shock using Taylor-McColl's eqn'ss approx
        %solution.
        
        R = 287.035;       
        %theta_c1= half angle of the nose cone. The rest of the symbols have their usual meanings
        
        
        %Approx value of shock angle using approx Taylor-McColl's Equation
        beta_c1 = (theta_c1 + (theta_c1^2 + ((154 * g_0^3 + 810 * g_0^2 + 990 * g_0 + 350)/(72 * (g_0+1)^4 * m_0^2)))^0.5)/((11 * g_0^2 + 50 * g_0 + 35)/(12 *(g_0 + 1)^2));

        %Wedge angle equivalent to conacal shock angle of beta_c1
        delta_1 = atand(2 * cotd(beta_c1) * (m_0^2 * (sind(beta_c1))^2 - 1)/(m_0^2 * (g_0 + cosd(2 * beta_c1)) + 2));

     mn1 = m_0* sind(beta_c1);
     
      %mach number after shock
     m_nfb = ((1 + ((g_0 - 1) * 0.5) * mn1^2)/(g_0 * mn1^2 - 0.2 * (g_0 - 1)))^0.5;
     m_cfb = m_nfb / sind(beta_c1 - delta_1);
     
     
     r_cfb = r_0 * (((g_0 + 1) * mn1^2)/(2 + (g_0 - 1) * mn1^2));
     p_cfb = p_0 * (1 + (2 * g_0 * (mn1^2 - 1))/(g_0 + 1));
     t_cfb = (p_cfb / p_0) * (r_0 / r_cfb) * t_0;
     v_cfb = m_cfb *(g_0 * R * t_cfb)^0.5;
     
     
     %return values
     co = [beta_c1, m_cfb, r_cfb, p_cfb, t_cfb, v_cfb];

end