function hy = eachReflection (th1, t_pre, v_pre, r_pre, h_tot, t_ref, a1) 

%function to calculate properties after each foreshock

%t_pre = temperature after previous shock
%v_pre = velocity after previous shock
%r_pre = density after previous shock

f = @(t) cpinteg(t,a1);
R = 287.035; %gas constant

            b_j = th1 : 0.1 : 89;
            z = length(b_j);
            it(1) = 1;
                    
            for j = 1:z
                             b = b_j(j);
                             s_b = ((R * t_pre + (sind(b) * v_pre)^2) * (tand(b - th1)/tand(b)) - (v_pre * cosd(b) * tand(b - th1))^2)/R;
                             
                             c = integral(f ,t_ref,s_b,'ArrayValued',true);
                             it(j+1) = 2 * (h_tot - c) * ((cosd(b - th1))^2) - (v_pre * cosd(b))^2;
                             
                            if it(j+1) * it(j) <0
                                
                                b_ob = b;
                                break;
                                
                            end
            end
                    
                 
            hy(1) = b_ob; %beta of the current shock
            
            %calculating properties if there wasnt unstart
            
            %temperature after shock reflection
            hy(2) = ((R * t_pre + (sind(b_ob) * v_pre)^2) * (tand(b_ob - th1)/tand(b_ob)) - (v_pre * cosd(b_ob) * tand(b_ob - th1))^2)/R; 
        
            
            hy(3) = v_pre * cosd(b_ob)/cosd(b_ob - th1);                 %velocity after each shock reflection
            hy(4) = r_pre * tand(b_ob)/tand(b_ob - th1);                 %density after shock 
                                        
            hy(5) = cpinteg(hy(2),a1);                                   %cp after each shock reflection
            hy(6) = hy(5)/(hy(5)- R);                                    %gamma after each shock reflection
                        
            hy(7) = hy(3)/(hy(6) * hy(2) * R)^0.5;                       %mach number after each shock reflection
            hy(8) = hy(4) * R * hy(2);                                   %pressure after each shock reflection  
            
end