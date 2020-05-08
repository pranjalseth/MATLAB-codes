function [tf, vf, rf, pf, bf, cpf, gf, mf]  = forebody(theta_n, aoa, t_0, v_0, r_0, h_tot, t_ref, e,a1)

    %function to return propertirs after the forebody
    %all the variabes mean the same as said in the main code
      
    R = 287.035;      % Gas constant
    f = @(t) cpinteg(t,a1); 
    
        th1 = (theta_n+aoa);
        
        t = t_0;
        v = v_0;
        r = r_0;
        
  for i=1:e        %controlling loop for the number of ramps  
        
        if e==1
            th1 = theta_n+aoa;
        else
            th1 = theta_n(i)-theta_n(i-1);
        end
        
        b_j = th1 : 0.1 : 89;
        z = length(b_j);
        it(1) = 1;
        
    for j = 1 : z  %loop to calculate the value of beta
        
           b = b_j(j);
           s_b = ((R * t + (sind(b) * v)^2)*(tand(b - th1) / tand(b)) - (v * cosd(b) * tand(b - th1))^2) /R;
           c = integral(f ,t_ref,s_b,'ArrayValued',true);
           it(j+1) = 2 * (h_tot - c) * ((cosd(b - th1))^2) - (v_0 * cosd(b))^2;
           
                            if it(j+1) * it(j)<0
                                
                                b_fb(i)=b;
                                break;
                                
                            end
    end
  
           
            t_fb(i) = ((R * t + (sind(b_fb(i)) * v)^2)*(tand(b_fb(i) - th1)/tand(b_fb(i))) - (v * cosd(b_fb) * tand(b_fb - th1))^2)/R;  %after foreshock
            r_fb(i) = r*tand(b_fb(i)) / tand(b_fb(i) - th1);     %density after foreshock
            
            v_fb(i) = v*cosd(b_fb(i)) / cosd(b_fb(i) - th1);     %velocity after foreshock
            p_fb(i) = r_fb(i) * R * t_fb(i);                           %pressure after foreshock
            
            t = t_fb(i);
            v = v_fb(i);
            p = p_fb(i);
            
            cp_fb(i) = cpinteg(t_fb(i), a1);
            g_fb(i) = cp_fb(i) /(cp_fb(i) - R);           %gamma after weach foreeshock 
            m_fb(i) = v_fb(i)/(g_fb(i) * t_fb(i) * R)^0.5;   %mach number after each foreshock 
  end
  
  tf = t_fb;
  vf = v_fb;
  rf = r_fb;
  pf = p_fb;
  bf = b_fb;
  cpf = cp_fb;
  gf = g_fb;
  mf = m_fb;
end

