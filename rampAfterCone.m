function [tf, vf, rf, pf, bf, mf] = rampAfterCone(g_0, e,theta_n, t_cfb, p_cfb, r_cfb, m_cfb)

      %function to return propertirs after the forebody
    %all the variabes mean the same as said in the main code
    
    %tf = temp matrix after foreshock(s) 
    %vf = velocity matrix after foreshock(s)
    %rf = density matrix after foreshock(s)
    %bf = pressure matrix after foreshock(s)
    %bf = shock angle(beta) matrix after foreshock(s)
    %gf = gamma matrix after foreshock(s)
    %mf = mach matrix after foreshock(s)
    %cpf = cp matrix after foreshock(s)
      
    R = 287.035;      % Gas constant

    
      
        
    t = t_cfb;
    r = r_cfb;
    m = m_cfb;
    p = p_cfb;
        
 %controlling loop for the number of ramps  
        
        it(1) = 1;
        
  
  for i = 1 : e     %loop to control the ramps
              
      th1 = theta_n(i+1)-theta_n(i);
      b_j = th1 : 0.1 : 89;
      z = length(b_j);
      
        b_j=th1:0.1:89;
        z=length(b_j);
        
                    for j=1:z
                    
                            a = 2*cotd(b_j(j)) * ((m * sind(b_j(j)))^2-1)/(m^2*(g_0+cosd(2*b_j(j)))+2);
                            b1 = tand(th1);
                            if (abs(a-b1)<= 0.001)
                            
                                b0=b_j(j);  
                                break;
                                
                            end
                    end
                    
                    b_fb(i) = b0;

    
                    mn1 = m*sind(b0);
                    m_n2 = ((1 + ((g_0-1)*0.5)*mn1^2)/(g_0*mn1^2-0.2*(g_0-1)))^0.5;

                    r_fb(i) = r * ((g_0+1)*mn1^2)/(2+(g_0-1)*mn1^2);

                    p_fb(i) = p *(1+(2*g_0/(g_0+1))*(mn1^2-1));

                    t_fb(i) = t * (p_fb(i)/p)*(r/r_fb(i));

                    m_fb(i) = m_n2 / sind(b0-th1);

                    v_fb(i) = m_fb(i) / (g_0 * t_fb(i) * R)^0.5;   %mach number after each foreshock


                    
                    t = t_fb(i);
                    r = r_fb(i);
                    m = m_fb(i);
                    p = p_fb(i);
  end
  
  
tf = t_fb;
vf = v_fb;
rf = r_fb;
pf = p_fb;
bf = b_fb;
mf = m_fb;




end

